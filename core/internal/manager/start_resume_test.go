package manager

import (
	"context"
	"encoding/json"
	"errors"
	"net/http"
	"os"
	"path/filepath"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager/resume"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func localSettings(t *testing.T, runDir, runID string) config.Settings {
	t.Helper()
	cfg := config.Defaults()
	cfg.Mode = config.ModeLocal
	cfg.RunID = runID
	cfg.RunDir = runDir
	cfg.RecordBatch = 16
	cfg.RecordInterval = 0.01
	return cfg
}

type mockCloudState struct {
	mu          sync.Mutex
	columnKeys  []string
	scalarSteps []int64
	saveNames   []string
	houseDown   atomic.Bool
	scalarCh    chan int64
}

func setupResumeCloud(t *testing.T) (*clienttest.Cloud, *mockCloudState) {
	t.Helper()
	cloud := clienttest.NewCloud(t)
	st := &mockCloudState{scalarCh: make(chan int64, 64)}

	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 1)
	cloud.ServeJSON(http.MethodGet, "/api/health", http.StatusOK, map[string]any{"ok": true})
	cloud.Handle(http.MethodGet, "/api/house/health", func(w http.ResponseWriter, _ *http.Request) {
		if st.houseDown.Load() {
			clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]any{"ok": false})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-1")
	cloud.ServeJSON(http.MethodPut, "/api/project/alice/demo/runs/exp-1/state", http.StatusOK, map[string]any{"ok": true})

	cloud.Handle(http.MethodPost, "/api/projects/alice/demo/series", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			Series []struct {
				Key string `json:"key"`
			} `json:"series"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		st.mu.Lock()
		for _, item := range body.Series {
			st.columnKeys = append(st.columnKeys, item.Key)
		}
		st.mu.Unlock()
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, r *http.Request) {
		if st.houseDown.Load() {
			clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]any{"error": "down"})
			return
		}
		var body struct {
			Type    string            `json:"type"`
			Metrics []json.RawMessage `json:"metrics"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		if body.Type == "scalar" {
			for _, raw := range body.Metrics {
				var entry struct {
					Index int64 `json:"index"`
				}
				if err := json.Unmarshal(raw, &entry); err == nil {
					st.mu.Lock()
					st.scalarSteps = append(st.scalarSteps, entry.Index)
					st.mu.Unlock()
					select {
					case st.scalarCh <- entry.Index:
					default:
					}
				}
			}
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	cloud.Handle(http.MethodPost, "/api/experiment/exp-1/files/prepare", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			Files []struct {
				Path string `json:"path"`
			} `json:"files"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		urls := make([]string, 0, len(body.Files))
		st.mu.Lock()
		for _, f := range body.Files {
			st.saveNames = append(st.saveNames, f.Path)
			urls = append(urls, cloud.URL()+"/upload")
		}
		st.mu.Unlock()
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"urls": urls})
	})
	cloud.ServeJSON(http.MethodPut, "/upload", http.StatusOK, map[string]any{"ok": true})
	cloud.ServeJSON(http.MethodPost, "/api/experiment/exp-1/files/complete", http.StatusCreated, map[string]any{"ok": true})

	return cloud, st
}

func waitForScalarCount(t *testing.T, st *mockCloudState, want int) {
	t.Helper()
	deadline := time.NewTimer(5 * time.Second)
	defer deadline.Stop()
	for {
		st.mu.Lock()
		n := len(st.scalarSteps)
		st.mu.Unlock()
		if n >= want {
			return
		}
		select {
		case <-st.scalarCh:
		case <-deadline.C:
			t.Fatalf("scalar count = %d, want >= %d", n, want)
		}
	}
}

func writeTestRecords(t *testing.T, path string, recs ...*recordv1.Record) int64 {
	t.Helper()
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	w, err := store.Create(path)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	for _, rec := range recs {
		payload, mErr := proto.Marshal(rec)
		if mErr != nil {
			_ = w.Close()
			t.Fatalf("Marshal: %v", mErr)
		}
		if err = w.Write(payload); err != nil {
			_ = w.Close()
			t.Fatalf("Write: %v", err)
		}
	}
	endOff := w.Offset()
	if err = w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return endOff
}

func startReq(runID string) *runv1.StartRecord {
	return &runv1.StartRecord{Project: "demo", Workspace: "alice", Id: runID}
}

func makeStartRecord(runID, workspace, project string) *recordv1.Record {
	return &recordv1.Record{
		Num: record.NumStart,
		RecordType: &recordv1.Record_Start{
			Start: &runv1.StartRecord{
				Id:        runID,
				Workspace: workspace,
				Project:   project,
			},
		},
	}
}

func scalarInput(step int64, val float64) *datav1.ScalarRecord {
	return &datav1.ScalarRecord{
		Key:       "loss",
		Step:      step,
		Type:      columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		Value:     &datav1.ScalarValue{Number: val},
		Timestamp: timestamppb.Now(),
	}
}

func upsertScalars(t *testing.T, mgr *Manager, sessionID string, recs ...*datav1.ScalarRecord) {
	t.Helper()
	if err := mgr.UpsertScalars(context.Background(), sessionID, recs); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
}

func makeColumnRecord(num int64, key string) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Column{
			Column: &columnv1.ColumnRecord{
				ColumnKey:  key,
				ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
			},
		},
	}
}

func makeScalarRecord(num int64, key string, step int64, val float64) *recordv1.Record {
	sc := scalarInput(step, val)
	sc.Key = key
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Scalar{Scalar: sc},
	}
}

func makeEndSaveRecord(num int64, name, sourcePath string) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{
				Name:       name,
				SourcePath: sourcePath,
				Policy:     savev1.SavePolicy_SAVE_POLICY_END,
				Type:       savev1.SaveType_SAVE_TYPE_CUSTOM,
			},
		},
	}
}

func TestResumeLocalModeAppendAndFinishGuard(t *testing.T) {
	runDir := filepath.Join(t.TempDir(), "run")
	cfg := localSettings(t, runDir, "local-resume")
	ctx := context.Background()

	// 第一次运行：写入 3 个 scalar (step 0..2) 与 1 条日志后“崩溃”退出（不写 Finish）。
	mgr1 := newTestManager()
	res1, err := mgr1.DeliverRunStart(ctx, cfg, startReq(cfg.RunID))
	if err != nil {
		t.Fatalf("first DeliverRunStart: %v", err)
	}
	if !res1.NewExperiment || res1.GlobalStep != -1 {
		t.Fatalf("first StartResult = %+v, want NewExperiment=true GlobalStep=-1", res1)
	}
	upsertScalars(t, mgr1, res1.SessionID, scalarInput(0, 1.0), scalarInput(1, 0.8), scalarInput(2, 0.6))
	if err = mgr1.UpsertLogs(ctx, res1.SessionID, []*terminalv1.LogRecord{{Line: "hello"}}); err != nil {
		t.Fatalf("first UpsertLogs: %v", err)
	}
	s1, err := mgr1.session(res1.SessionID)
	if err != nil {
		t.Fatalf("session1: %v", err)
	}
	prevCounter := s1.counter.Value()
	prevEpoch := s1.epoch.Value()
	if err = mgr1.Shutdown(ctx); err != nil {
		t.Fatalf("mgr1.Shutdown: %v", err)
	}

	// 同目录重启：重建指标、计数器与 epoch，不重复写 Start。
	mgr2 := newTestManager()
	res2, err := mgr2.DeliverRunStart(ctx, cfg, startReq(cfg.RunID))
	if err != nil {
		t.Fatalf("second DeliverRunStart: %v", err)
	}
	if res2.NewExperiment {
		t.Fatal("resumed local StartResult.NewExperiment = true, want false")
	}
	if res2.GlobalStep != 2 {
		t.Fatalf("resumed GlobalStep = %d, want 2", res2.GlobalStep)
	}
	s2, err := mgr2.session(res2.SessionID)
	if err != nil {
		t.Fatalf("session2: %v", err)
	}
	if got := s2.counter.Value(); got != prevCounter {
		t.Fatalf("resumed counter = %d, want %d", got, prevCounter)
	}
	if got := s2.epoch.Value(); got != prevEpoch {
		t.Fatalf("resumed epoch = %d, want %d", got, prevEpoch)
	}

	// 重复 step=2 被跳过；step=3 正常写入且不再重建 loss 列。
	upsertScalars(t, mgr2, res2.SessionID, scalarInput(2, 9.9), scalarInput(3, 0.4))
	if err = mgr2.DeliverRunFinish(ctx, res2.SessionID, &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err = mgr2.Confirm(ctx, res2.SessionID); err != nil {
		t.Fatalf("Confirm: %v", err)
	}

	recs := readStore(t, runFilePath(t, cfg))
	var starts, finishes, cols, scalars int
	for _, r := range recs {
		switch record.TypeOf(r) {
		case record.TypeStart:
			starts++
		case record.TypeFinish:
			finishes++
		case record.TypeColumn:
			cols++
		case record.TypeScalar:
			scalars++
		}
	}
	if starts != 1 || finishes != 1 || cols != 1 || scalars != 4 {
		t.Fatalf("record counts = (start=%d finish=%d col=%d scalar=%d), want (1, 1, 1, 4)",
			starts, finishes, cols, scalars)
	}

	// 已含 FinishRecord 的 run 文件拒绝再次续传。
	mgr3 := newTestManager()
	_, err = mgr3.DeliverRunStart(ctx, cfg, startReq(cfg.RunID))
	if !errors.Is(err, ErrInvalidStartRecord) {
		t.Fatalf("third start after finish err = %v, want ErrInvalidStartRecord", err)
	}
}

func TestResumeTruncatesCorruptTailAndPreservesExistingOnRollback(t *testing.T) {
	runDir := filepath.Join(t.TempDir(), "run")
	cfg := localSettings(t, runDir, "corrupt-tail")
	ctx := context.Background()

	mgr1 := newTestManager()
	res1, err := mgr1.DeliverRunStart(ctx, cfg, startReq(cfg.RunID))
	if err != nil {
		t.Fatalf("first DeliverRunStart: %v", err)
	}
	upsertScalars(t, mgr1, res1.SessionID, scalarInput(0, 1.0), scalarInput(1, 0.5))
	if err = mgr1.Shutdown(ctx); err != nil {
		t.Fatalf("mgr1.Shutdown: %v", err)
	}

	path := runFilePath(t, cfg)
	// 在文件末尾附加 5 字节垃圾，模拟崩溃时的半条记录。
	fp, err := os.OpenFile(path, os.O_WRONLY|os.O_APPEND, 0o644)
	if err != nil {
		t.Fatalf("OpenFile: %v", err)
	}
	if _, err = fp.Write([]byte{0xDE, 0xAD, 0xBE, 0xEF, 0x42}); err != nil {
		_ = fp.Close()
		t.Fatalf("append garbage: %v", err)
	}
	_ = fp.Close()

	// 续传回滚时不得删除已存在的 run 文件。
	rbRun, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	rebuilt, err := resume.LoadLocalRebuild(cfg, rbRun)
	if err != nil {
		t.Fatalf("loadLocalRebuild: %v", err)
	}
	sRollback := newSession(sessionOptions{})
	sRollback.ctx, sRollback.cancel = context.WithCancel(ctx)
	sRollback.cfg = cfg
	sRollback.run = rbRun
	if err = sRollback.open(sessionInit{cfg: cfg, run: rbRun, resume: rebuilt}); err != nil {
		t.Fatalf("sRollback.open: %v", err)
	}
	sRollback.rollback()
	if _, statErr := os.Stat(path); statErr != nil {
		t.Fatalf("rollback must preserve existing run file, got stat err=%v", statErr)
	}

	// 正常重启应截断尾部并继续写入。
	mgr2 := newTestManager()
	res2, err := mgr2.DeliverRunStart(ctx, cfg, startReq(cfg.RunID))
	if err != nil {
		t.Fatalf("second DeliverRunStart: %v", err)
	}
	if res2.GlobalStep != 1 {
		t.Fatalf("GlobalStep = %d, want 1", res2.GlobalStep)
	}
	upsertScalars(t, mgr2, res2.SessionID, scalarInput(2, 0.25))
	if err = mgr2.DeliverRunFinish(ctx, res2.SessionID, &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err = mgr2.Confirm(ctx, res2.SessionID); err != nil {
		t.Fatalf("Confirm: %v", err)
	}

	recs := readStore(t, path)
	var scalars int
	for _, r := range recs {
		if record.TypeOf(r) == record.TypeScalar {
			scalars++
		}
	}
	if scalars != 3 {
		t.Fatalf("scalars = %d, want 3 after corrupt tail truncation + append", scalars)
	}
}

func TestResumeOnlineCheckpointBackfillAndNoDuplicateColumns(t *testing.T) {
	cloud, st := setupResumeCloud(t)
	cfg := startSettings(t, cloud)
	cfg.RunID = "online-resume"
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	ctx := context.Background()

	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	runFile, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}
	cpFile, err := run.CheckpointFile()
	if err != nil {
		t.Fatalf("CheckpointFile: %v", err)
	}

	// 构造本地 store：seq 1..4 (Start + Col(loss) + Scalar(0..1)) 已进入 checkpoint，
	// seq 5..6 (Scalar(2..3)) 为崩溃前尚未确认的尾部记录。
	records := []*recordv1.Record{
		makeStartRecord(cfg.RunID, "alice", "demo"),
		makeColumnRecord(1, "loss"),
		makeScalarRecord(2, "loss", 0, 1.0),
		makeScalarRecord(3, "loss", 1, 0.8),
		makeScalarRecord(4, "loss", 2, 0.6),
		makeScalarRecord(5, "loss", 3, 0.4),
	}
	writeTestRecords(t, runFile, records...)
	cpOff, err := store.EndOffset(runFile, nil, 4, -1)
	if err != nil {
		t.Fatalf("EndOffset: %v", err)
	}
	if err = store.WriteCheckpoint(cpFile, store.Checkpoint{
		Version:   store.CheckpointVersion,
		RunID:     cfg.RunID,
		StoreSeq:  4,
		Offset:    cpOff,
		UpdatedAt: time.Now(),
	}); err != nil {
		t.Fatalf("WriteCheckpoint: %v", err)
	}

	mgr := newTestManager()
	res, err := mgr.DeliverRunStart(ctx, cfg, startReq(cfg.RunID))
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	if res.GlobalStep != 3 {
		t.Fatalf("GlobalStep = %d, want 3 from local rebuild", res.GlobalStep)
	}

	// 续传新提交 step=4（不再重建 loss 列），并等待 seq 5..6 与新记录全部送达。
	upsertScalars(t, mgr, res.SessionID, scalarInput(4, 0.2))
	waitForScalarCount(t, st, 3) // step 2, 3 (从 checkpoint 补投) + step 4 (新提交)

	if err = mgr.DeliverRunFinish(ctx, res.SessionID, &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err = mgr.Confirm(ctx, res.SessionID); err != nil {
		t.Fatalf("Confirm: %v", err)
	}

	st.mu.Lock()
	cols := append([]string(nil), st.columnKeys...)
	steps := append([]int64(nil), st.scalarSteps...)
	st.mu.Unlock()
	if len(cols) != 0 {
		t.Fatalf("expected 0 column uploads after checkpoint before seq=2, got %v", cols)
	}
	wantSteps := map[int64]bool{2: true, 3: true, 4: true}
	if len(steps) != len(wantSteps) {
		t.Fatalf("uploaded scalar steps = %v, want only [2 3 4]", steps)
	}
	for _, step := range steps {
		if !wantSteps[step] {
			t.Fatalf("unexpected uploaded step %d in %v", step, steps)
		}
	}

	// 确认最终 store 中 start 仍只有 1 条、finish 1 条、scalar 共 5 条 (step 0..4)。
	all := readStore(t, runFile)
	var starts, finishes, scalars int
	for _, r := range all {
		switch record.TypeOf(r) {
		case record.TypeStart:
			starts++
		case record.TypeFinish:
			finishes++
		case record.TypeScalar:
			scalars++
		}
	}
	if starts != 1 || finishes != 1 || scalars != 5 {
		t.Fatalf("store counts = (start=%d finish=%d scalar=%d), want (1, 1, 5)", starts, finishes, scalars)
	}
}

func TestResumeMissingAndOutOfBoundsCheckpointFallback(t *testing.T) {
	cases := []struct {
		name  string
		write func(t *testing.T, cpFile, runID string, endOff int64)
	}{
		{
			name:  "missing checkpoint",
			write: func(_ *testing.T, _ string, _ string, _ int64) {},
		},
		{
			name: "seq beyond records",
			write: func(t *testing.T, cpFile, runID string, endOff int64) {
				if err := store.WriteCheckpoint(cpFile, store.Checkpoint{
					Version:  store.CheckpointVersion,
					RunID:    runID,
					StoreSeq: 99,
					Offset:   endOff,
				}); err != nil {
					t.Fatal(err)
				}
			},
		},
		{
			name: "offset beyond end",
			write: func(t *testing.T, cpFile, runID string, endOff int64) {
				if err := store.WriteCheckpoint(cpFile, store.Checkpoint{
					Version:  store.CheckpointVersion,
					RunID:    runID,
					StoreSeq: 2,
					Offset:   endOff + 100,
				}); err != nil {
					t.Fatal(err)
				}
			},
		},
		{
			name: "misaligned boundary offset",
			write: func(t *testing.T, cpFile, runID string, _ int64) {
				if err := store.WriteCheckpoint(cpFile, store.Checkpoint{
					Version:  store.CheckpointVersion,
					RunID:    runID,
					StoreSeq: 2,
					Offset:   store.HeaderLen + 2,
				}); err != nil {
					t.Fatal(err)
				}
			},
		},
		{
			name:  "valid boundary belonging to another store seq",
			write: func(_ *testing.T, _ string, _ string, _ int64) {},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud, st := setupResumeCloud(t)
			cfg := startSettings(t, cloud)
			cfg.RunID = "cp-fallback"
			cfg.RunDir = filepath.Join(t.TempDir(), "run")
			ctx := context.Background()

			run, err := runctx.New(cfg)
			if err != nil {
				t.Fatalf("runctx.New: %v", err)
			}
			runFile, err := run.RunFile()
			if err != nil {
				t.Fatalf("RunFile: %v", err)
			}
			cpFile, err := run.CheckpointFile()
			if err != nil {
				t.Fatalf("CheckpointFile: %v", err)
			}
			endOff := writeTestRecords(t, runFile,
				makeStartRecord(cfg.RunID, "alice", "demo"),
				makeColumnRecord(1, "loss"),
				makeScalarRecord(2, "loss", 0, 1.0),
				makeScalarRecord(3, "loss", 1, 0.5),
			)
			if tc.name == "valid boundary belonging to another store seq" {
				wrongOff, endErr := store.EndOffset(runFile, nil, 3, -1)
				if endErr != nil {
					t.Fatalf("EndOffset: %v", endErr)
				}
				if wErr := store.WriteCheckpoint(cpFile, store.Checkpoint{
					Version:  store.CheckpointVersion,
					RunID:    cfg.RunID,
					StoreSeq: 2,
					Offset:   wrongOff,
				}); wErr != nil {
					t.Fatal(wErr)
				}
			} else {
				tc.write(t, cpFile, cfg.RunID, endOff)
			}

			mgr := newTestManager()
			res, err := mgr.DeliverRunStart(ctx, cfg, &runv1.StartRecord{
				Project:   "demo",
				Workspace: "alice",
				Id:        cfg.RunID,
			})
			if err != nil {
				t.Fatalf("DeliverRunStart: %v", err)
			}
			waitForScalarCount(t, st, 2)

			if err = mgr.DeliverRunFinish(ctx, res.SessionID, &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
				t.Fatalf("DeliverRunFinish: %v", err)
			}
			if _, err = mgr.Confirm(ctx, res.SessionID); err != nil {
				t.Fatalf("Confirm: %v", err)
			}

			st.mu.Lock()
			cols := append([]string(nil), st.columnKeys...)
			steps := append([]int64(nil), st.scalarSteps...)
			st.mu.Unlock()
			if len(cols) != 1 || cols[0] != "loss" {
				t.Fatalf("fallback columnKeys = %v, want [loss]", cols)
			}
			if len(steps) != 2 {
				t.Fatalf("fallback scalarSteps = %v, want 2 records from head", steps)
			}
		})
	}
}

func TestResumeOnlineMergesCloudSummaryAndEndSaves(t *testing.T) {
	cloud, st := setupResumeCloud(t)
	// 返回已存在的实验与云端摘要，验证云端高水位与 EndSaves 合并恢复。
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusOK, map[string]any{
		"cuid":      "exp-1",
		"slug":      "run-1",
		"name":      "existing-exp",
		"createdAt": "2026-01-01T00:00:00Z",
	})
	cloud.ServeJSON(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-1", http.StatusOK, map[string]any{
		"scalar": []map[string]any{
			{"key": "loss", "step": 10},
			{"key": "acc", "step": 5},
		},
		"log": []map[string]any{{"step": 20}},
	})

	cfg := startSettings(t, cloud)
	cfg.RunID = "merge-resume"
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	endFile := filepath.Join(t.TempDir(), "end.txt")
	if err := os.WriteFile(endFile, []byte("final-artifact"), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	ctx := context.Background()

	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	runFile, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}
	cpFile, err := run.CheckpointFile()
	if err != nil {
		t.Fatalf("CheckpointFile: %v", err)
	}

	// 本地只有 step=4 的 loss 与一条 END 保存记录，全部进入 checkpoint。
	endOff := writeTestRecords(t, runFile,
		makeStartRecord(cfg.RunID, "alice", "demo"),
		makeColumnRecord(1, "loss"),
		makeScalarRecord(2, "loss", 4, 0.9),
		makeEndSaveRecord(3, "end.txt", endFile),
	)
	if err = store.WriteCheckpoint(cpFile, store.Checkpoint{
		Version:  store.CheckpointVersion,
		RunID:    cfg.RunID,
		StoreSeq: 4,
		Offset:   endOff,
	}); err != nil {
		t.Fatalf("WriteCheckpoint: %v", err)
	}

	mgr := newTestManager()
	res, err := mgr.DeliverRunStart(ctx, cfg, &runv1.StartRecord{
		Project:   "demo",
		Workspace: "alice",
		Id:        cfg.RunID,
		Resume:    runv1.ResumeMode_RESUME_MODE_ALLOW,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	if res.GlobalStep != 10 {
		t.Fatalf("merged GlobalStep = %d, want 10 (cloud higher)", res.GlobalStep)
	}
	s, err := mgr.session(res.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	if got := s.epoch.Value(); got != 20 {
		t.Fatalf("merged epoch = %d, want 20", got)
	}
	if pre := s.metrics.PrecheckScalar(&datav1.ScalarRecord{
		Key:  "loss",
		Step: 8,
		Type: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
	}); pre.NeedsColumn || pre.Decision != metrics.DecisionSkip {
		t.Fatalf("loss step 8 precheck = %+v, want NeedsColumn=false Skip", pre)
	}

	// 新写 step=11 并结束：此前本地的 END 保存记录也必须在 Close 时上传。
	upsertScalars(t, mgr, res.SessionID, scalarInput(11, 0.1))
	waitForScalarCount(t, st, 1)

	if err = mgr.DeliverRunFinish(ctx, res.SessionID, &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err = mgr.Confirm(ctx, res.SessionID); err != nil {
		t.Fatalf("Confirm: %v", err)
	}

	st.mu.Lock()
	saves := append([]string(nil), st.saveNames...)
	st.mu.Unlock()
	if len(saves) != 1 || saves[0] != "end.txt" {
		t.Fatalf("uploaded END saves = %v, want [end.txt]", saves)
	}
}

func TestResumeSkipStoreIgnoresExistingRunFile(t *testing.T) {
	cloud, _ := setupResumeCloud(t)
	cfg := startSettings(t, cloud)
	cfg.RunID = "skip-resume"
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	cfg.SkipStore = true

	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	runFile, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}
	// 即使 RunDir 里放了一个已结束的 run 文件，skip_store 也不读不续传。
	writeTestRecords(t, runFile, backfillStart(), makeColumnRecord(1, "loss"), makeScalarRecord(2, "loss", 5, 1.0), backfillFinish())

	mgr := newTestManager()
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Project:   "demo",
		Workspace: "alice",
		Id:        cfg.RunID,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart(skip_store): %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })
	if res.GlobalStep != -1 {
		t.Fatalf("skip_store GlobalStep = %d, want -1", res.GlobalStep)
	}
}

func TestResumeOnlineUploadsAllHistoryAfterLocalCrash(t *testing.T) {
	runDir := filepath.Join(t.TempDir(), "run")
	localCfg := localSettings(t, runDir, "local-to-online")
	ctx := context.Background()

	// 1. local 模式写入历史记录后崩溃退出（不写 Finish），过程中不得产生 checkpoint.json。
	mgr1 := newTestManager()
	res1, err := mgr1.DeliverRunStart(ctx, localCfg, startReq(localCfg.RunID))
	if err != nil {
		t.Fatalf("local DeliverRunStart: %v", err)
	}
	upsertScalars(t, mgr1, res1.SessionID, scalarInput(0, 1.0), scalarInput(1, 0.5))
	if err = mgr1.Shutdown(ctx); err != nil {
		t.Fatalf("local Shutdown: %v", err)
	}

	run, err := runctx.New(localCfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	cpFile, err := run.CheckpointFile()
	if err != nil {
		t.Fatalf("CheckpointFile: %v", err)
	}
	if _, statErr := os.Stat(cpFile); !os.IsNotExist(statErr) {
		t.Fatalf("local run must not create checkpoint file %s, got stat err=%v", cpFile, statErr)
	}

	// 2. 同目录切换 online 续传：必须从头回放并上传全部 local 历史记录，而非跳过。
	cloud, st := setupResumeCloud(t)
	onlineCfg := startSettings(t, cloud)
	onlineCfg.RunID = localCfg.RunID
	onlineCfg.RunDir = runDir

	mgr2 := newTestManager()
	res2, err := mgr2.DeliverRunStart(ctx, onlineCfg, startReq(onlineCfg.RunID))
	if err != nil {
		t.Fatalf("online DeliverRunStart: %v", err)
	}
	upsertScalars(t, mgr2, res2.SessionID, scalarInput(2, 0.25))
	waitForScalarCount(t, st, 3)

	if err = mgr2.DeliverRunFinish(ctx, res2.SessionID, &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err = mgr2.Confirm(ctx, res2.SessionID); err != nil {
		t.Fatalf("Confirm: %v", err)
	}

	st.mu.Lock()
	cols := append([]string(nil), st.columnKeys...)
	steps := append([]int64(nil), st.scalarSteps...)
	st.mu.Unlock()
	if len(cols) != 1 || cols[0] != "loss" {
		t.Fatalf("uploaded columns = %v, want [loss]", cols)
	}
	wantSteps := map[int64]bool{0: true, 1: true, 2: true}
	if len(steps) != len(wantSteps) {
		t.Fatalf("uploaded scalar steps = %v, want all local history + new [0 1 2]", steps)
	}
	for _, step := range steps {
		if !wantSteps[step] {
			t.Fatalf("unexpected uploaded step %d in %v", step, steps)
		}
	}
}

func TestResumeRejectsInvalidOrMismatchedStartRecord(t *testing.T) {
	cases := []struct {
		name string
		recs []*recordv1.Record
		req  *runv1.StartRecord
	}{
		{
			name: "missing start record in non-empty store",
			recs: []*recordv1.Record{makeColumnRecord(1, "loss"), makeScalarRecord(2, "loss", 0, 1.0)},
			req:  startReq("resume-check"),
		},
		{
			name: "start record not at head",
			recs: []*recordv1.Record{makeScalarRecord(1, "loss", 0, 1.0), makeStartRecord("resume-check", "alice", "demo")},
			req:  startReq("resume-check"),
		},
		{
			name: "start record with wrong num",
			recs: []*recordv1.Record{{
				Num:        1,
				RecordType: &recordv1.Record_Start{Start: startReq("resume-check")},
			}},
			req: startReq("resume-check"),
		},
		{
			name: "run id mismatch",
			recs: []*recordv1.Record{makeStartRecord("other-run", "alice", "demo"), makeScalarRecord(1, "loss", 0, 1.0)},
			req:  startReq("resume-check"),
		},
		{
			name: "project mismatch",
			recs: []*recordv1.Record{makeStartRecord("resume-check", "alice", "other-proj"), makeScalarRecord(1, "loss", 0, 1.0)},
			req:  startReq("resume-check"),
		},
		{
			name: "workspace mismatch after normalization",
			recs: []*recordv1.Record{makeStartRecord("resume-check", "bob", "demo"), makeScalarRecord(1, "loss", 0, 1.0)},
			req:  startReq("resume-check"),
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud, st := setupResumeCloud(t)
			cfg := startSettings(t, cloud)
			cfg.RunID = "resume-check"
			cfg.RunDir = filepath.Join(t.TempDir(), "run")

			run, err := runctx.New(cfg)
			if err != nil {
				t.Fatalf("runctx.New: %v", err)
			}
			runFile, err := run.RunFile()
			if err != nil {
				t.Fatalf("RunFile: %v", err)
			}
			endOff := writeTestRecords(t, runFile, tc.recs...)

			mgr := newTestManager()
			if _, startErr := mgr.DeliverRunStart(context.Background(), cfg, tc.req); !errors.Is(startErr, ErrInvalidStartRecord) {
				t.Fatalf("DeliverRunStart err = %v, want ErrInvalidStartRecord", startErr)
			}
			// 校验失败时不得创建实验、修改本地文件或触发回放上传。
			if createReqs := cloud.RequestsFor(http.MethodPost, "/api/project/alice/demo/experiment"); len(createReqs) != 0 {
				t.Fatalf("CreateOrResumeExperiment requests = %d, want 0", len(createReqs))
			}
			info, statErr := os.Stat(runFile)
			if statErr != nil || info.Size() != endOff {
				t.Fatalf("runFile modified: statErr=%v size=%d want %d", statErr, info.Size(), endOff)
			}
			st.mu.Lock()
			uploaded := len(st.scalarSteps)
			st.mu.Unlock()
			if uploaded != 0 {
				t.Fatalf("uploaded scalars = %d, want 0", uploaded)
			}
		})
	}
}
