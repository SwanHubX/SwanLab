package resume

import (
	"errors"
	"os"
	"path/filepath"
	"testing"

	"go.uber.org/goleak"
	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

func writeRecords(t *testing.T, path string, recs ...*recordv1.Record) int64 {
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

func startRec() *recordv1.Record {
	return &recordv1.Record{Num: record.NumStart, RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}}}
}

func scalarRec(num int64, step int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Scalar{
			Scalar: &datav1.ScalarRecord{
				Key:   "loss",
				Step:  step,
				Type:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
				Value: &datav1.ScalarValue{Number: float64(step)},
			},
		},
	}
}

func finishRec() *recordv1.Record {
	return &recordv1.Record{Num: record.NumFinish, RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{}}}
}

func endSaveRec(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{
				Name:   "end.txt",
				Policy: savev1.SavePolicy_SAVE_POLICY_END,
				Type:   savev1.SaveType_SAVE_TYPE_CUSTOM,
			},
		},
	}
}

func TestLoadLocalRebuildCases(t *testing.T) {
	cfg := config.Defaults()
	cfg.RunID = "run-1"
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}

	// 不存在时返回 (nil, nil)
	res, err := LoadLocalRebuild(cfg, run)
	if err != nil || res != nil {
		t.Fatalf("missing run file = (%v, %v), want (nil, nil)", res, err)
	}

	runFile, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}

	// 已结束运行拒绝续传
	writeRecords(t, runFile, startRec(), scalarRec(1, 0), finishRec())
	if _, err = LoadLocalRebuild(cfg, run); !errors.Is(err, ErrInvalidRunFile) {
		t.Fatalf("finished run err = %v, want ErrInvalidRunFile", err)
	}

	// 缺失 StartRecord、Start 不在首部、Start Num 非法均拒绝续传
	for _, tc := range []struct {
		name string
		recs []*recordv1.Record
	}{
		{name: "missing start", recs: []*recordv1.Record{scalarRec(1, 0)}},
		{name: "start not first", recs: []*recordv1.Record{scalarRec(1, 0), startRec()}},
		{name: "wrong start num", recs: []*recordv1.Record{{Num: 1, RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}}}}},
	} {
		_ = os.Remove(runFile)
		writeRecords(t, runFile, tc.recs...)
		if _, err = LoadLocalRebuild(cfg, run); !errors.Is(err, ErrInvalidRunFile) {
			t.Fatalf("%s err = %v, want ErrInvalidRunFile", tc.name, err)
		}
	}

	// 只有文件头的空文件允许继续（视作尚未写入 StartRecord）
	_ = os.Remove(runFile)
	writeRecords(t, runFile)
	if res, err = LoadLocalRebuild(cfg, run); err != nil || res == nil || res.Records != 0 {
		t.Fatalf("empty header-only file = (%+v, %v), want 0 records", res, err)
	}

	// SkipStore 忽略本地文件
	skipCfg := cfg
	skipCfg.SkipStore = true
	if res, err = LoadLocalRebuild(skipCfg, run); err != nil || res != nil {
		t.Fatalf("skip store = (%v, %v), want (nil, nil)", res, err)
	}

	// 未结束运行可成功重建
	_ = os.Remove(runFile)
	writeRecords(t, runFile, startRec(), scalarRec(1, 0))
	res, err = LoadLocalRebuild(cfg, run)
	if err != nil || res == nil || res.Records != 2 {
		t.Fatalf("valid rebuild = (%+v, %v), want 2 records", res, err)
	}
}

func TestValidateStartIdentity(t *testing.T) {
	base := &metrics.RebuildResult{
		Records:  2,
		HasStart: true,
		Start:    &runv1.StartRecord{Id: "run-1", Workspace: "alice", Project: "demo"},
	}

	if err := ValidateStart(nil, StartIdentity{RunID: "run-1", Project: "demo"}); err != nil {
		t.Fatalf("nil rebuild = %v, want nil", err)
	}
	if err := ValidateStart(&metrics.RebuildResult{Records: 0}, StartIdentity{RunID: "run-1", Project: "demo"}); err != nil {
		t.Fatalf("empty rebuild = %v, want nil", err)
	}
	if err := ValidateStart(base, StartIdentity{RunID: "run-1", Workspace: "alice", Project: "demo"}); err != nil {
		t.Fatalf("matching start = %v, want nil", err)
	}

	for _, tc := range []struct {
		name   string
		local  *metrics.RebuildResult
		target StartIdentity
	}{
		{
			name:   "run id mismatch",
			local:  base,
			target: StartIdentity{RunID: "run-2", Workspace: "alice", Project: "demo"},
		},
		{
			name:   "project mismatch",
			local:  base,
			target: StartIdentity{RunID: "run-1", Workspace: "alice", Project: "other"},
		},
		{
			name:   "workspace mismatch",
			local:  base,
			target: StartIdentity{RunID: "run-1", Workspace: "bob", Project: "demo"},
		},
		{
			name: "local empty workspace mismatch when strict",
			local: &metrics.RebuildResult{
				Records:  2,
				HasStart: true,
				Start:    &runv1.StartRecord{Id: "run-1", Project: "demo"},
			},
			target: StartIdentity{RunID: "run-1", Workspace: "alice", Project: "demo"},
		},
	} {
		if err := ValidateStart(tc.local, tc.target); !errors.Is(err, ErrInvalidRunFile) {
			t.Fatalf("%s err = %v, want ErrInvalidRunFile", tc.name, err)
		}
	}

	// AllowEmptyWorkspace 允许早期 local(workspace="") 在转在线时被规范化补齐。
	localNoSpace := &metrics.RebuildResult{
		Records:  2,
		HasStart: true,
		Start:    &runv1.StartRecord{Id: "run-1", Project: "demo"},
	}
	if err := ValidateStart(localNoSpace, StartIdentity{
		RunID: "run-1", Workspace: "alice", Project: "demo", AllowEmptyWorkspace: true,
	}); err != nil {
		t.Fatalf("local empty workspace into online err = %v, want nil", err)
	}
}

func TestResolveCheckpointAndCountUploadable(t *testing.T) {
	cfg := config.Defaults()
	cfg.RunID = "run-1"
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	runFile, _ := run.RunFile()
	cpFile, _ := run.CheckpointFile()

	endOff := writeRecords(t, runFile, startRec(), scalarRec(1, 0), endSaveRec(2), scalarRec(3, 1))
	local, err := metrics.Rebuild(runFile)
	if err != nil {
		t.Fatalf("Rebuild: %v", err)
	}

	// 缺失 checkpoint 回退文件头
	seq, off, valid := ResolveCheckpoint(cpFile, run.RunID, runFile, local)
	if valid || seq != 0 || off != store.HeaderLen {
		t.Fatalf("missing cp = (%d, %d, %v), want (0, %d, false)", seq, off, valid, store.HeaderLen)
	}

	// 统计整文件可立即上传记录数：Start/END_SAVE 跳过，2 条 scalar 计入
	if got := CountUploadable(runFile, store.HeaderLen, endOff); got != 2 {
		t.Fatalf("CountUploadable = %d, want 2", got)
	}

	// 合法中段 checkpoint：前 2 条记录结束 offset
	cpOff, err := store.EndOffset(runFile, nil, 2, -1)
	if err != nil {
		t.Fatalf("EndOffset: %v", err)
	}
	if err = store.WriteCheckpoint(cpFile, store.Checkpoint{
		Version:  store.CheckpointVersion,
		RunID:    run.RunID,
		StoreSeq: 2,
		Offset:   cpOff,
	}); err != nil {
		t.Fatalf("WriteCheckpoint: %v", err)
	}

	seq, off, valid = ResolveCheckpoint(cpFile, run.RunID, runFile, local)
	if !valid || seq != 2 || off != cpOff {
		t.Fatalf("valid cp = (%d, %d, %v), want (2, %d, true)", seq, off, valid, cpOff)
	}

	// 非记录边界 offset 应回退
	_ = store.WriteCheckpoint(cpFile, store.Checkpoint{
		Version:  store.CheckpointVersion,
		RunID:    run.RunID,
		StoreSeq: 2,
		Offset:   store.HeaderLen + 2,
	})
	seq, off, valid = ResolveCheckpoint(cpFile, run.RunID, runFile, local)
	if valid || seq != 0 || off != store.HeaderLen {
		t.Fatalf("bad boundary cp = (%d, %d, %v), want fallback", seq, off, valid)
	}

	// 合法记录边界但属于其它 store_seq 应拒绝并回退
	wrongSeqOff, err := store.EndOffset(runFile, nil, 3, -1)
	if err != nil {
		t.Fatalf("EndOffset(3): %v", err)
	}
	_ = store.WriteCheckpoint(cpFile, store.Checkpoint{
		Version:  store.CheckpointVersion,
		RunID:    run.RunID,
		StoreSeq: 2,
		Offset:   wrongSeqOff,
	})
	seq, off, valid = ResolveCheckpoint(cpFile, run.RunID, runFile, local)
	if valid || seq != 0 || off != store.HeaderLen {
		t.Fatalf("wrong-seq boundary cp = (%d, %d, %v), want fallback", seq, off, valid)
	}

	// InitOnlineProgress 根据有效 checkpoint 设定重放起点
	_ = store.WriteCheckpoint(cpFile, store.Checkpoint{
		Version:  store.CheckpointVersion,
		RunID:    run.RunID,
		StoreSeq: 2,
		Offset:   cpOff,
	})
	prog := InitOnlineProgress(ProgressOptions{
		Mode:           config.ModeOnline,
		RunID:          run.RunID,
		RunFile:        runFile,
		CheckpointFile: cpFile,
		HasStore:       true,
		Resume:         local,
	})
	if !prog.CheckpointValid || prog.CheckpointSeq != 2 || prog.ReplayFrom != 3 || prog.ReplayStartOffset != cpOff {
		t.Fatalf("InitOnlineProgress = %+v, want seq=2 replayFrom=3 off=%d", prog, cpOff)
	}
}

func TestOpenStoreCreateAndAppend(t *testing.T) {
	cfg := config.Defaults()
	cfg.RunID = "run-open"
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}

	opened, err := OpenStore(StoreOptions{Run: run, Queue: 16})
	if err != nil {
		t.Fatalf("OpenStore(create): %v", err)
	}
	if !opened.Created || opened.Writer == nil || opened.RunFile == "" || opened.CheckpointFile == "" {
		t.Fatalf("unexpected OpenedStore: %+v", opened)
	}
	_ = opened.Writer.Close()
	_ = os.Remove(opened.RunFile)

	writeRecords(t, opened.RunFile, startRec(), scalarRec(1, 0))
	local, err := metrics.Rebuild(opened.RunFile)
	if err != nil {
		t.Fatalf("Rebuild: %v", err)
	}
	appended, err := OpenStore(StoreOptions{Run: run, Resume: local, Queue: 16})
	if err != nil {
		t.Fatalf("OpenStore(append): %v", err)
	}
	if appended.Created {
		t.Fatal("append OpenedStore.Created = true, want false")
	}
	_ = appended.Writer.Close()

	skipped, err := OpenStore(StoreOptions{SkipStore: true, Run: run, Queue: 16})
	if err != nil || skipped.Writer != nil || skipped.Created {
		t.Fatalf("SkipStore OpenedStore = %+v, err=%v", skipped, err)
	}
}
