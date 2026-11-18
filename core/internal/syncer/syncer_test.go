package syncer

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

	"github.com/google/uuid"
	"go.uber.org/goleak"
	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

type recordedCloud struct {
	cloud *clienttest.Cloud

	mu         sync.Mutex
	expReqs    []map[string]any
	columns    []map[string]any
	metrics    []map[string]any
	logs       []map[string]any
	profiles   []map[string]any
	stopStates []string
	summaryReq atomic.Int32
}

func newRecordedCloud(t *testing.T, newExperiment bool, summary map[string]any) *recordedCloud {
	t.Helper()
	cloud := clienttest.NewCloud(t)
	rc := &recordedCloud{cloud: cloud}

	cloud.ServeJSON(http.MethodGet, "/api/auth/verify", http.StatusOK, map[string]any{
		"uid":       1,
		"username":  "alice",
		"createdAt": "2026-01-01T00:00:00Z",
	})
	for _, proj := range []struct{ owner, name string }{
		{"alice", "demo"},
		{"orig-ws", "orig-proj"},
	} {
		owner, name := proj.owner, proj.name
		cloud.ServeJSON(http.MethodGet, "/api/project/"+owner+"/"+name, http.StatusOK, map[string]any{
			"cuid":       "proj-1",
			"name":       name,
			"path":       "/@alice/" + name,
			"visibility": "PRIVATE",
			"version":    1,
			"_count":     map[string]any{"experiments": 3},
			"group":      map[string]any{"username": "alice"},
		})
	}
	cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		rc.mu.Lock()
		rc.expReqs = append(rc.expReqs, body)
		rc.mu.Unlock()

		status := http.StatusOK
		if newExperiment {
			status = http.StatusCreated
		}
		clienttest.WriteJSON(w, status, map[string]any{
			"cuid":      "exp-1",
			"slug":      "exp-slug",
			"name":      body["name"],
			"createdAt": "2026-01-01T00:00:00Z",
		})
	})
	cloud.Handle(http.MethodPost, "/api/project/alice/orig-proj/experiment", func(w http.ResponseWriter, r *http.Request) {
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		rc.mu.Lock()
		rc.expReqs = append(rc.expReqs, body)
		rc.mu.Unlock()

		status := http.StatusOK
		if newExperiment {
			status = http.StatusCreated
		}
		clienttest.WriteJSON(w, status, map[string]any{
			"cuid":      "exp-1",
			"slug":      "exp-slug",
			"name":      body["name"],
			"createdAt": "2026-01-01T00:00:00Z",
		})
	})
	for _, proj := range []string{"demo", "orig-proj"} {
		cloud.Handle(http.MethodPut, "/api/project/alice/"+proj+"/runs/exp-1/state", func(w http.ResponseWriter, r *http.Request) {
			var body map[string]any
			_ = json.NewDecoder(r.Body).Decode(&body)
			rc.mu.Lock()
			if state, ok := body["state"].(string); ok {
				rc.stopStates = append(rc.stopStates, state)
			}
			rc.mu.Unlock()
			clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
		})
		cloud.Handle(http.MethodPut, "/api/project/alice/"+proj+"/runs/exp-1/profile", func(w http.ResponseWriter, r *http.Request) {
			var body map[string]any
			_ = json.NewDecoder(r.Body).Decode(&body)
			rc.mu.Lock()
			rc.profiles = append(rc.profiles, body)
			rc.mu.Unlock()
			clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
		})
		cloud.Handle(http.MethodPost, "/api/projects/alice/"+proj+"/series", func(w http.ResponseWriter, r *http.Request) {
			var body map[string]any
			_ = json.NewDecoder(r.Body).Decode(&body)
			rc.mu.Lock()
			if cols, ok := body["series"].([]any); ok {
				for _, c := range cols {
					if m, ok := c.(map[string]any); ok {
						rc.columns = append(rc.columns, m)
					}
				}
			}
			rc.mu.Unlock()
			clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
		})
	}
	cloud.Handle(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-1", func(w http.ResponseWriter, _ *http.Request) {
		rc.summaryReq.Add(1)
		if summary == nil {
			summary = map[string]any{}
		}
		clienttest.WriteJSON(w, http.StatusOK, summary)
	})
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, r *http.Request) {
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		rc.mu.Lock()
		if body["type"] == "log" {
			rc.logs = append(rc.logs, body)
		} else {
			rc.metrics = append(rc.metrics, body)
		}
		rc.mu.Unlock()
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cloud.Handle(http.MethodPost, "/api/resources/presigned/put", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			Paths []string `json:"paths"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		urls := make([]string, 0, len(body.Paths))
		for _, p := range body.Paths {
			urls = append(urls, cloud.URL()+"/upload/"+p)
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"urls": urls})
	})
	cloud.Handle(http.MethodPut, "/upload/media/image/pic.png", func(w http.ResponseWriter, _ *http.Request) {
		w.WriteHeader(http.StatusOK)
	})
	return rc
}

func baseSyncSettings(runDir string, cloud *clienttest.Cloud) config.Settings {
	cfg := config.Defaults()
	cfg.RunDir = runDir
	cfg.RunID = "target-run"
	cfg.RecordBatch = 16
	cfg.RecordInterval = 0.01
	cfg.APIKey = "test-key"
	if cloud != nil {
		cfg.APIHost = cloud.URL()
		cfg.WebHost = cloud.URL()
	}
	return cfg
}

func writeSyncRunFile(t *testing.T, runDir, filename string, records ...*recordv1.Record) string {
	t.Helper()
	if err := os.MkdirAll(runDir, 0o755); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(runDir, filename)
	w, err := store.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, rec := range records {
		b, mErr := proto.Marshal(rec)
		if mErr != nil {
			t.Fatal(mErr)
		}
		if err := w.Write(b); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

func startRecord(id string) *recordv1.Record {
	return &recordv1.Record{
		Num: record.NumStart,
		RecordType: &recordv1.Record_Start{
			Start: &runv1.StartRecord{
				Id:        id,
				Workspace: "orig-ws",
				Project:   "orig-proj",
				Name:      "orig-name",
				Color:     "#123456",
			},
		},
	}
}

func scalarMetricRecord(num int64, key string, step int64, val float64, at *timestamppb.Timestamp) *recordv1.Record {
	return &recordv1.Record{
		Num:       num,
		Timestamp: at,
		RecordType: &recordv1.Record_Scalar{
			Scalar: &datav1.ScalarRecord{
				Key:       key,
				Step:      step,
				Type:      columnv1.ColumnType_COLUMN_TYPE_SCALAR,
				Value:     &datav1.ScalarValue{Number: val},
				Timestamp: at,
			},
		},
	}
}

func waitSyncDrained(t *testing.T, mgr *Manager, sid string) {
	t.Helper()
	deadline := time.Now().Add(5 * time.Second)
	ticker := time.NewTicker(10 * time.Millisecond)
	defer ticker.Stop()
	for {
		st, err := mgr.Stats(sid)
		if err != nil {
			t.Fatalf("Stats: %v", err)
		}
		if st.GetState() == operationv1.CoreState_CORE_STATE_FINISHED {
			return
		}
		if time.Now().After(deadline) {
			t.Fatalf("timed out waiting for FINISHED stats, last=%v", st)
		}
		<-ticker.C
	}
}

func TestDeliverSyncStartValidation(t *testing.T) {
	mgr := NewManager(Options{})
	defer func() { _ = mgr.Close() }()

	// 1. 缺 run_dir
	if _, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: config.Defaults()}); !errors.Is(err, ErrInvalidSettings) {
		t.Fatalf("missing run_dir got %v, want ErrInvalidSettings", err)
	}

	// 2. 空目录（无 run-*.swanlab）
	emptyDir := t.TempDir()
	if _, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(emptyDir, nil)}); !errors.Is(err, ErrInvalidStartRecord) {
		t.Fatalf("missing run file got %v, want ErrInvalidStartRecord", err)
	}

	// 3. 多个 run-*.swanlab
	multiDir := t.TempDir()
	writeSyncRunFile(t, multiDir, "run-1.swanlab", startRecord("1"))
	writeSyncRunFile(t, multiDir, "run-2.swanlab", startRecord("2"))
	if _, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(multiDir, nil)}); !errors.Is(err, ErrInvalidStartRecord) {
		t.Fatalf("multiple run files got %v, want ErrInvalidStartRecord", err)
	}

	// 4. 首条非 StartRecord
	badStartDir := t.TempDir()
	writeSyncRunFile(t, badStartDir, "run-bad.swanlab", scalarMetricRecord(1, "loss", 1, 0.5, ts(100, 0)))
	if _, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(badStartDir, nil)}); !errors.Is(err, ErrInvalidStartRecord) {
		t.Fatalf("non-start first record got %v, want ErrInvalidStartRecord", err)
	}

	// 5. 非法 SyncRange
	validDir := t.TempDir()
	writeSyncRunFile(t, validDir, "run-ok.swanlab", startRecord("ok"))
	if _, err := mgr.DeliverSyncStart(context.Background(), StartOptions{
		Settings: baseSyncSettings(validDir, nil),
		Range:    &corev1.SyncRange{StartStep: i64(10), EndStep: i64(5)},
	}); !errors.Is(err, ErrInvalidSyncRange) {
		t.Fatalf("invalid range got %v, want ErrInvalidSyncRange", err)
	}

	// 6. 离线环境变量 / skip_store=true 下仍允许启动
	cfg := baseSyncSettings(validDir, nil)
	cfg.Mode = config.ModeOffline
	cfg.SkipStore = true
	cfg.RunID = "" // 未传 run_id 时从 startRecord 回填
	sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: cfg})
	if err != nil {
		t.Fatalf("DeliverSyncStart with offline/skip_store: %v", err)
	}
	if err := mgr.Disconnect(context.Background(), sid); err != nil {
		t.Fatalf("Disconnect: %v", err)
	}
}

func TestSyncFullLifecycleEndToEnd(t *testing.T) {
	rc := newRecordedCloud(t, true, nil)
	runDir := t.TempDir()
	mediaDir := filepath.Join(runDir, "media", "image")
	if err := os.MkdirAll(mediaDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(mediaDir, "pic.png"), []byte("pngbytes"), 0o644); err != nil {
		t.Fatal(err)
	}
	filesDir := filepath.Join(runDir, "files")
	if err := os.MkdirAll(filesDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(filesDir, "config.yaml"), []byte("lr: 0.01\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	xAxis := "epoch"
	hidden := true
	writeSyncRunFile(t, runDir, "run-orig-1.swanlab",
		startRecord("orig-1"),
		&recordv1.Record{
			Num: 1,
			RecordType: &recordv1.Record_Column{
				Column: &columnv1.ColumnRecord{
					ColumnKey:   "train/loss",
					ColumnType:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
					SectionName: "Training",
					XAxis:       &xAxis,
					Hidden:      &hidden,
				},
			},
		},
		scalarMetricRecord(2, "train/loss", 1, 0.5, ts(100, 0)),
		scalarMetricRecord(3, "train/loss", 1, 0.4, ts(101, 0)), // 重复 step，应跳过
		scalarMetricRecord(4, "acc", 2, 0.9, ts(102, 0)),        // 自动建列
		&recordv1.Record{
			Num: 5,
			RecordType: &recordv1.Record_Media{
				Media: &datav1.MediaRecord{
					Key:       "img",
					Step:      1,
					Type:      columnv1.ColumnType_COLUMN_TYPE_IMAGE,
					Timestamp: ts(103, 0),
					Value: &datav1.MediaValue{
						Items: []*datav1.MediaItem{{Filename: "pic.png", Caption: "sample"}},
					},
				},
			},
		},
		&recordv1.Record{
			Num: 6,
			RecordType: &recordv1.Record_Histogram{
				Histogram: &datav1.HistogramRecord{
					Key:       "weights",
					Step:      1,
					Type:      columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM,
					Timestamp: ts(104, 0),
					Value:     &datav1.HistogramValue{Values: []uint64{1, 2}, Bins: []float64{0, 1, 2}},
				},
			},
		},
		&recordv1.Record{
			Num: 7,
			RecordType: &recordv1.Record_Log{
				Log: &terminalv1.LogRecord{
					Line:      "hello world",
					Level:     terminalv1.LogLevel_LOG_LEVEL_INFO,
					Timestamp: ts(105, 0),
					Epoch:     1,
				},
			},
		},
		&recordv1.Record{
			Num: record.NumConfigSave,
			RecordType: &recordv1.Record_Save{
				Save: &savev1.SaveRecord{
					Name:       "config.yaml",
					SourcePath: filepath.Join(filesDir, "config.yaml"),
					Type:       savev1.SaveType_SAVE_TYPE_CONFIG,
				},
			},
		},
		&recordv1.Record{
			Num: record.NumFinish,
			RecordType: &recordv1.Record_Finish{
				Finish: &runv1.FinishRecord{
					State:      runv1.RunState_RUN_STATE_FINISHED,
					FinishedAt: ts(110, 0),
				},
			},
		},
	)

	mgr := NewManager(Options{})
	defer func() { _ = mgr.Close() }()

	sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{
		Settings:  baseSyncSettings(runDir, rc.cloud),
		Workspace: "alice",
		Project:   "demo",
		RunID:     "override-id",
	})
	if err != nil {
		t.Fatalf("DeliverSyncStart: %v", err)
	}
	// DeliverSyncStart 不应发起云端请求。
	rc.mu.Lock()
	if len(rc.expReqs) != 0 {
		rc.mu.Unlock()
		t.Fatalf("DeliverSyncStart made unexpected cloud calls: %+v", rc.expReqs)
	}
	rc.mu.Unlock()

	flushRes, err := mgr.DeliverSyncFlush(context.Background(), sid)
	if err != nil {
		t.Fatalf("DeliverSyncFlush: %v", err)
	}
	if flushRes.Path != "/@alice/demo/exp-slug" {
		t.Fatalf("Path = %q, want /@alice/demo/exp-slug", flushRes.Path)
	}

	waitSyncDrained(t, mgr, sid)
	confirmRes, err := mgr.ConfirmSyncFinish(context.Background(), sid)
	if err != nil {
		t.Fatalf("ConfirmSyncFinish: %v", err)
	}
	if !confirmRes.Success || confirmRes.Message != confirmOKMessage {
		t.Fatalf("ConfirmSyncFinish = %+v, want Success=true", confirmRes)
	}

	rc.mu.Lock()
	defer rc.mu.Unlock()
	if len(rc.expReqs) != 1 || rc.expReqs[0]["cuid"] != "override-id" {
		t.Fatalf("expReqs = %+v, want cuid=override-id", rc.expReqs)
	}
	if rc.summaryReq.Load() != 0 {
		t.Fatalf("new experiment should not fetch summary, got %d", rc.summaryReq.Load())
	}
	if len(rc.columns) != 4 {
		t.Fatalf("uploaded columns = %d, want 4 (%+v)", len(rc.columns), rc.columns)
	}
	if len(rc.metrics) != 3 {
		t.Fatalf("metric batches = %d, want 3 (scalar + media + histogram)", len(rc.metrics))
	}
	if len(rc.logs) != 1 {
		t.Fatalf("log batches = %d, want 1", len(rc.logs))
	}
	if len(rc.profiles) != 1 {
		t.Fatalf("profiles = %d, want 1", len(rc.profiles))
	}
	if len(rc.stopStates) != 1 || rc.stopStates[0] != "FINISHED" {
		t.Fatalf("stopStates = %v, want [FINISHED]", rc.stopStates)
	}
}

func TestSyncIncrementalExistingExperimentAndMissingFinishRecord(t *testing.T) {
	rc := newRecordedCloud(t, false, map[string]any{
		"log":       []map[string]any{{"key": "stdout", "step": 2}},
		"scalar":    []map[string]any{{"key": "loss", "step": 2}},
		"histogram": []map[string]any{{"key": "weights", "step": 3}},
	})
	runDir := t.TempDir()
	writeSyncRunFile(t, runDir, "run-inc.swanlab",
		startRecord("inc"),
		scalarMetricRecord(1, "loss", 1, 0.8, ts(100, 0)), // <= summary(2) 跳过
		scalarMetricRecord(2, "loss", 3, 0.2, ts(101, 0)), // 上传
		&recordv1.Record{
			Num: 3,
			RecordType: &recordv1.Record_Histogram{
				Histogram: &datav1.HistogramRecord{
					Key: "weights", Step: 3, Type: columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM,
					Timestamp: ts(102, 0),
					Value:     &datav1.HistogramValue{Values: []uint64{1}, Bins: []float64{0, 1}},
				},
			},
		}, // <= summary(3) 跳过
		&recordv1.Record{
			Num: 4,
			RecordType: &recordv1.Record_Log{
				Log: &terminalv1.LogRecord{Line: "old", Epoch: 2, Level: terminalv1.LogLevel_LOG_LEVEL_INFO, Timestamp: ts(103, 0)},
			},
		}, // <= epoch(2) 跳过
		&recordv1.Record{
			Num: 5,
			RecordType: &recordv1.Record_Log{
				Log: &terminalv1.LogRecord{Line: "new", Epoch: 3, Level: terminalv1.LogLevel_LOG_LEVEL_INFO, Timestamp: ts(104, 0)},
			},
		}, // 上传
	)

	mgr := NewManager(Options{})
	defer func() { _ = mgr.Close() }()

	sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{
		Settings: baseSyncSettings(runDir, rc.cloud),
	})
	if err != nil {
		t.Fatal(err)
	}
	if _, flushErr := mgr.DeliverSyncFlush(context.Background(), sid); flushErr != nil {
		t.Fatal(flushErr)
	}
	waitSyncDrained(t, mgr, sid)
	res, err := mgr.ConfirmSyncFinish(context.Background(), sid)
	if err != nil {
		t.Fatal(err)
	}
	if !res.Success {
		t.Fatalf("ConfirmSyncFinish = %+v", res)
	}

	rc.mu.Lock()
	defer rc.mu.Unlock()
	if rc.summaryReq.Load() != 1 {
		t.Fatalf("summaryReq = %d, want 1", rc.summaryReq.Load())
	}
	// 摘要已定义 loss 与 weights，不应重复建列。
	if len(rc.columns) != 0 {
		t.Fatalf("columns = %+v, want 0", rc.columns)
	}
	if len(rc.metrics) != 1 {
		t.Fatalf("metric batches = %d, want 1 (only loss step 3)", len(rc.metrics))
	}
	if len(rc.logs) != 1 {
		t.Fatalf("log batches = %d, want 1", len(rc.logs))
	}
	// 已有实验且无 FinishRecord：保持 RUNNING，不上报 StopExperiment。
	if len(rc.stopStates) != 0 {
		t.Fatalf("stopStates = %v, want empty", rc.stopStates)
	}
}

func TestSyncNewExperimentWithoutFinishMarksCrashedAndAbortedUploadsErrorLog(t *testing.T) {
	t.Run("new without finish marks crashed without fake error log", func(t *testing.T) {
		rc := newRecordedCloud(t, true, nil)
		runDir := t.TempDir()
		writeSyncRunFile(t, runDir, "run-nofin.swanlab",
			startRecord("nofin"),
			scalarMetricRecord(1, "loss", 1, 0.5, ts(100, 0)),
		)
		mgr := NewManager(Options{})
		defer func() { _ = mgr.Close() }()

		sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(runDir, rc.cloud)})
		if err != nil {
			t.Fatal(err)
		}
		if _, flushErr := mgr.DeliverSyncFlush(context.Background(), sid); flushErr != nil {
			t.Fatal(flushErr)
		}
		waitSyncDrained(t, mgr, sid)
		res, err := mgr.ConfirmSyncFinish(context.Background(), sid)
		if err != nil || !res.Success {
			t.Fatalf("confirm = (%+v, %v)", res, err)
		}

		rc.mu.Lock()
		defer rc.mu.Unlock()
		if len(rc.logs) != 0 {
			t.Fatalf("unexpected logs: %+v", rc.logs)
		}
		if len(rc.stopStates) != 1 || rc.stopStates[0] != "CRASHED" {
			t.Fatalf("stopStates = %v, want [CRASHED]", rc.stopStates)
		}
	})

	t.Run("failed finish record uploads error log and reports state", func(t *testing.T) {
		rc := newRecordedCloud(t, false, nil)
		runDir := t.TempDir()
		writeSyncRunFile(t, runDir, "run-abort.swanlab",
			startRecord("abort"),
			scalarMetricRecord(1, "loss", 1, 0.5, ts(100, 0)),
			&recordv1.Record{
				Num: record.NumFinish,
				RecordType: &recordv1.Record_Finish{
					Finish: &runv1.FinishRecord{
						State:      runv1.RunState_RUN_STATE_ABORTED,
						Error:      "user aborted",
						FinishedAt: ts(110, 0),
					},
				},
			},
		)
		mgr := NewManager(Options{})
		defer func() { _ = mgr.Close() }()

		sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(runDir, rc.cloud)})
		if err != nil {
			t.Fatal(err)
		}
		if _, flushErr := mgr.DeliverSyncFlush(context.Background(), sid); flushErr != nil {
			t.Fatal(flushErr)
		}
		waitSyncDrained(t, mgr, sid)
		res, err := mgr.ConfirmSyncFinish(context.Background(), sid)
		if err != nil || !res.Success {
			t.Fatalf("confirm = (%+v, %v)", res, err)
		}

		rc.mu.Lock()
		defer rc.mu.Unlock()
		if len(rc.logs) != 1 {
			t.Fatalf("logs = %d, want 1", len(rc.logs))
		}
		if len(rc.stopStates) != 1 || rc.stopStates[0] != "ABORTED" {
			t.Fatalf("stopStates = %v, want [ABORTED]", rc.stopStates)
		}
	})
}

func TestSyncRangeReplayAllowsOldStepsAndFiltersBeforeAutoColumn(t *testing.T) {
	rc := newRecordedCloud(t, false, map[string]any{
		"scalar": []map[string]any{{"key": "loss", "step": 10}},
	})
	runDir := t.TempDir()
	writeSyncRunFile(t, runDir, "run-range.swanlab",
		startRecord("range"),
		scalarMetricRecord(1, "loss", 2, 0.5, ts(100, 0)),
		scalarMetricRecord(2, "loss", 5, 0.2, ts(200, 0)),
		scalarMetricRecord(3, "filtered_out", 99, 1.0, ts(300, 0)),
		&recordv1.Record{
			Num: 4,
			RecordType: &recordv1.Record_Log{
				Log: &terminalv1.LogRecord{Line: "skip-without-time-range", Epoch: 1, Level: terminalv1.LogLevel_LOG_LEVEL_INFO, Timestamp: ts(100, 0)},
			},
		},
		&recordv1.Record{
			Num: record.NumFinish,
			RecordType: &recordv1.Record_Finish{
				Finish: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_CRASHED, Error: "boom", FinishedAt: ts(300, 0)},
			},
		},
	)

	mgr := NewManager(Options{})
	defer func() { _ = mgr.Close() }()

	sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{
		Settings: baseSyncSettings(runDir, rc.cloud),
		Range:    &corev1.SyncRange{StartStep: i64(1), EndStep: i64(5)},
	})
	if err != nil {
		t.Fatal(err)
	}
	if _, flushErr := mgr.DeliverSyncFlush(context.Background(), sid); flushErr != nil {
		t.Fatal(flushErr)
	}
	waitSyncDrained(t, mgr, sid)
	res, err := mgr.ConfirmSyncFinish(context.Background(), sid)
	if err != nil || !res.Success {
		t.Fatalf("confirm = (%+v, %v)", res, err)
	}

	rc.mu.Lock()
	defer rc.mu.Unlock()
	// filtered_out 在过滤后被丢弃，不应自动建列。
	if len(rc.columns) != 0 {
		t.Fatalf("columns = %+v, want 0", rc.columns)
	}
	// step 2 & 5 <= summary(10)，但 AllowReplay 放行并打包进 1 个标量批次。
	if len(rc.metrics) != 1 {
		t.Fatalf("metric batches = %d, want 1", len(rc.metrics))
	}
	if rows, _ := rc.metrics[0]["metrics"].([]any); len(rows) != 2 {
		t.Fatalf("uploaded rows = %+v, want 2 rows", rc.metrics[0]["metrics"])
	}
	// 仅设置了 step 范围，终端日志与崩溃日志均不回放。
	if len(rc.logs) != 0 {
		t.Fatalf("logs = %+v, want 0", rc.logs)
	}
}

func TestSyncCorruptTailAndInvalidPayload(t *testing.T) {
	rc := newRecordedCloud(t, true, nil)
	runDir := t.TempDir()
	path := filepath.Join(runDir, "run-corrupt.swanlab")
	w, err := store.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, rec := range []*recordv1.Record{
		startRecord("corrupt"),
		scalarMetricRecord(1, "loss", 1, 0.5, ts(100, 0)),
	} {
		b, _ := proto.Marshal(rec)
		if wErr := w.Write(b); wErr != nil {
			t.Fatal(wErr)
		}
	}
	// 写入一条非法 protobuf 载荷（应跳过并继续）
	if wErr := w.Write([]byte{0x80}); wErr != nil {
		t.Fatal(wErr)
	}
	b2, _ := proto.Marshal(scalarMetricRecord(2, "loss", 2, 0.3, ts(101, 0)))
	if wErr := w.Write(b2); wErr != nil {
		t.Fatal(wErr)
	}
	if err = w.Close(); err != nil {
		t.Fatal(err)
	}

	// 追加半截坏 chunk 尾部（应保留前缀并停止扫描）
	fp, err := os.OpenFile(path, os.O_WRONLY|os.O_APPEND, 0o644)
	if err != nil {
		t.Fatal(err)
	}
	if _, wErr := fp.Write([]byte{0x01, 0x02, 0x03}); wErr != nil {
		t.Fatal(wErr)
	}
	_ = fp.Close()

	mgr := NewManager(Options{})
	defer func() { _ = mgr.Close() }()

	sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(runDir, rc.cloud)})
	if err != nil {
		t.Fatal(err)
	}
	if _, flushErr := mgr.DeliverSyncFlush(context.Background(), sid); flushErr != nil {
		t.Fatal(flushErr)
	}
	waitSyncDrained(t, mgr, sid)
	res, err := mgr.ConfirmSyncFinish(context.Background(), sid)
	if err != nil || !res.Success {
		t.Fatalf("confirm = (%+v, %v)", res, err)
	}

	rc.mu.Lock()
	defer rc.mu.Unlock()
	if len(rc.metrics) != 1 {
		t.Fatalf("metric batches = %d, want 1", len(rc.metrics))
	}
	if rows, _ := rc.metrics[0]["metrics"].([]any); len(rows) != 2 {
		t.Fatalf("uploaded rows = %+v, want 2 valid scalar rows", rc.metrics[0]["metrics"])
	}
}

func TestSyncDisconnectAndStopFailure(t *testing.T) {
	t.Run("disconnect does not call stop_experiment", func(t *testing.T) {
		rc := newRecordedCloud(t, true, nil)
		runDir := t.TempDir()
		writeSyncRunFile(t, runDir, "run-disc.swanlab",
			startRecord("disc"),
			scalarMetricRecord(1, "loss", 1, 0.5, ts(100, 0)),
		)
		mgr := NewManager(Options{})
		defer func() { _ = mgr.Close() }()

		sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(runDir, rc.cloud)})
		if err != nil {
			t.Fatal(err)
		}
		if _, flushErr := mgr.DeliverSyncFlush(context.Background(), sid); flushErr != nil {
			t.Fatal(flushErr)
		}
		if err := mgr.Disconnect(context.Background(), sid); err != nil {
			t.Fatalf("Disconnect: %v", err)
		}
		if _, err := mgr.Stats(sid); !errors.Is(err, ErrSessionNotFound) {
			t.Fatalf("Stats after Disconnect = %v, want ErrSessionNotFound", err)
		}
		rc.mu.Lock()
		stops := len(rc.stopStates)
		rc.mu.Unlock()
		if stops != 0 {
			t.Fatalf("Disconnect must not stop experiment, got %d stops", stops)
		}
	})

	t.Run("stop experiment error returns friendly message", func(t *testing.T) {
		rc := newRecordedCloud(t, true, nil)
		rc.cloud.ServeJSON(http.MethodPut, "/api/project/alice/orig-proj/runs/exp-1/state", http.StatusBadRequest, map[string]any{"error": "bad"})
		runDir := t.TempDir()
		writeSyncRunFile(t, runDir, "run-stoperr.swanlab",
			startRecord("stoperr"),
			scalarMetricRecord(1, "loss", 1, 0.5, ts(100, 0)),
			&recordv1.Record{
				Num:        record.NumFinish,
				RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}},
			},
		)
		mgr := NewManager(Options{})
		defer func() { _ = mgr.Close() }()

		sid, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(runDir, rc.cloud)})
		if err != nil {
			t.Fatal(err)
		}
		if _, flushErr := mgr.DeliverSyncFlush(context.Background(), sid); flushErr != nil {
			t.Fatal(flushErr)
		}
		waitSyncDrained(t, mgr, sid)
		res, err := mgr.ConfirmSyncFinish(context.Background(), sid)
		if err != nil {
			t.Fatalf("ConfirmSyncFinish err = %v", err)
		}
		if res.Success || res.Message != confirmFinishFailedMsg {
			t.Fatalf("ConfirmSyncFinish = %+v, want failure message %q", res, confirmFinishFailedMsg)
		}
	})
}

func TestSyncEnqueueOversizedRecordDoesNotBlock(t *testing.T) {
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	rec := scalarMetricRecord(1, "very_long_metric_key_exceeding_tiny_byte_budget", 1, 1.0, ts(100, 0))
	s := &syncSession{
		ctx:     ctx,
		budget:  budget.New(budget.Limits{Bytes: 4}),
		tracker: uploader.NewTracker(),
		inCh:    make(chan *recordv1.Record, 2),
	}
	if !s.enqueue(rec) {
		t.Fatal("enqueue should succeed for oversized single record via Adjust fallback")
	}
	usedRecs, usedBytes := s.budget.Usage()
	if usedRecs != 1 || usedBytes != record.Size(rec) {
		t.Fatalf("budget usage = (%d, %d), want (1, %d)", usedRecs, usedBytes, record.Size(rec))
	}
	if got := s.tracker.Snapshot().GetTotalRecords(); got != 1 {
		t.Fatalf("total_records = %d, want 1", got)
	}
}

func TestSyncSessionIDUsesTimeOrderedUUIDv7(t *testing.T) {
	runDir := t.TempDir()
	writeSyncRunFile(t, runDir, "run-uuidv7.swanlab", startRecord("u1"))

	mgr := NewManager(Options{})
	defer func() { _ = mgr.Close() }()

	first, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(runDir, nil)})
	if err != nil {
		t.Fatalf("first DeliverSyncStart: %v", err)
	}
	second, err := mgr.DeliverSyncStart(context.Background(), StartOptions{Settings: baseSyncSettings(runDir, nil)})
	if err != nil {
		t.Fatalf("second DeliverSyncStart: %v", err)
	}

	for _, raw := range []string{first, second} {
		parsed, pErr := uuid.Parse(raw)
		if pErr != nil {
			t.Fatalf("uuid.Parse(%q): %v", raw, pErr)
		}
		if parsed.Version() != 7 {
			t.Fatalf("uuid %q version = %d, want 7", raw, parsed.Version())
		}
	}
	if first >= second {
		t.Fatalf("expected UUIDv7 monotonicity: first=%q second=%q", first, second)
	}
}
