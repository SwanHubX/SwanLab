package manager

import (
	"context"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func newCheckpointTestSession(
	t *testing.T,
	mode config.Mode,
	sender uploader.Sender,
) (*Session, string, string) {
	t.Helper()
	cfg := config.Defaults()
	cfg.Mode = mode
	cfg.RunID = "cp-run"
	cfg.RunDir = t.TempDir()
	cfg.RecordBatch = 16
	cfg.RecordInterval = 0.005

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
	writer, err := store.Create(runFile)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{
		Queue:         64,
		FlushBytes:    1,
		FlushInterval: 2 * time.Millisecond,
		IndexStride:   1,
	})
	if err != nil {
		_ = writer.Close()
		t.Fatalf("NewAsyncWriter: %v", err)
	}

	b := budget.New(budget.Limits{Records: 64, Bytes: 1 << 20})
	tracker := uploader.NewTracker()
	var progress *progressTracker
	var onStored func(*recordv1.Record, uint64)
	if mode == config.ModeOnline {
		progress = newProgressTracker(0)
		onStored = progress.onStored
	}

	intake, err := buffer.New(buffer.Options{
		StoreOnly:   mode != config.ModeOnline,
		Store:       async.In(),
		UploadQueue: 64,
		Budget:      b,
		Flushed: func() uint64 {
			seq, _ := async.Flushed()
			return seq
		},
		OnAdmit:  tracker.AddTotalRecords,
		OnStored: onStored,
	})
	if err != nil {
		_ = async.Close()
		t.Fatalf("buffer.New: %v", err)
	}

	s := newSession(sessionOptions{})
	s.cfg = cfg
	s.run = run
	s.store = async
	s.intake = intake
	s.budget = b
	s.tracker = tracker
	s.metrics = metrics.New()
	s.counter = counter.New(0)
	s.epoch = counter.New(0)
	s.watcher = watcher.NullWatcher{}
	s.runFile = runFile
	s.runFileCreated = true
	s.checkpointFile = cpFile
	s.progress = progress

	if mode == config.ModeOnline {
		up, err := uploader.New(uploader.Options{
			In:     intake.UploadOut(),
			Budget: b,
			Senders: map[record.Type]uploader.Sender{
				record.TypeScalar: sender,
			},
			Batch:     16,
			Interval:  5 * time.Millisecond,
			Tracker:   tracker,
			OnConfirm: progress.resolve,
			OnDrop:    progress.resolveOne,
			OnDrained: func() {
				s.stopCheckpoint()
				s.flushCheckpoint()
				s.markDrained()
			},
		})
		if err != nil {
			_ = intake.Close(context.Background())
			_ = async.Close()
			t.Fatalf("uploader.New: %v", err)
		}
		s.uploader = up
		s.drainFn = up.Wait
	}
	return s, runFile, cpFile
}

func waitStoreFlushed(t *testing.T, s *store.AsyncWriter, want uint64) {
	t.Helper()
	deadline := time.Now().Add(2 * time.Second)
	for {
		seq, _ := s.Flushed()
		if seq >= want {
			return
		}
		if time.Now().After(deadline) {
			t.Fatalf("flushed seq = %d, want >= %d", seq, want)
		}
		runtimeGosched()
	}
}

func runtimeGosched() {
	timer := time.NewTimer(time.Millisecond)
	<-timer.C
}

func TestCheckpointBoundsAndFinalFlushOnline(t *testing.T) {
	gate := make(chan struct{})
	entered := make(chan struct{})
	collect := &collectingSender{notify: make(chan struct{}, 8)}
	gated := &gatedSender{inner: collect, gate: gate, entered: entered}

	s, runFile, cpFile := newCheckpointTestSession(t, config.ModeOnline, gated)
	t.Cleanup(func() {
		s.stopCheckpoint()
		if s.uploader != nil {
			_ = s.uploader.Close()
		}
		if s.store != nil {
			_ = s.store.Close()
		}
	})

	ticked := make(chan struct{}, 16)
	s.checkpointInterval = 5 * time.Millisecond
	s.checkpointTicked = func() {
		select {
		case ticked <- struct{}{}:
		default:
		}
	}
	s.startCheckpoint()

	// seq=1: Start (立即完结); seq=2: Scalar (将被 gatedSender 阻塞上传)。
	startRec := backfillStart()
	scRec := backfillScalar(1)
	commitToIntake(t, s.intake, s.cfg, startRec, scRec)

	waitStoreFlushed(t, s.store, 2)
	<-entered // uploader 已拿到 seq=2 并阻塞在 gate，此时 watermark=1 < flushedSeq=2。

	select {
	case <-ticked:
	case <-time.After(2 * time.Second):
		t.Fatal("timeout waiting for initial checkpoint tick")
	}

	cp, ok, err := store.LoadCheckpoint(cpFile, s.run.RunID)
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint(blocked) = (%v, %v), want (nil, true)", err, ok)
	}
	if cp.StoreSeq != 1 {
		t.Fatalf("StoreSeq while seq=2 unconfirmed = %d, want 1", cp.StoreSeq)
	}
	flushedSeq, flushedOff := s.store.Flushed()
	if flushedSeq != 2 {
		t.Fatalf("flushedSeq = %d, want 2", flushedSeq)
	}
	wantOff, err := store.EndOffset(runFile, s.store.Index(), 1, flushedOff)
	if err != nil {
		t.Fatalf("EndOffset: %v", err)
	}
	if cp.Offset != wantOff {
		t.Fatalf("Offset = %d, want %d", cp.Offset, wantOff)
	}

	// 放行上传并触发收尾：OnDrained 应停止协程并落盘最终连续水位（含 seq=3 Finish）。
	close(gate)
	collect.waitFor(t, 1)

	if err = s.finish(context.Background(), &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("finish: %v", err)
	}
	if err = s.waitDrained(context.Background()); err != nil {
		t.Fatalf("waitDrained: %v", err)
	}

	cpFinal, ok, err := store.LoadCheckpoint(cpFile, s.run.RunID)
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint(final) = (%v, %v), want (nil, true)", err, ok)
	}
	finalSeq, finalOff := s.store.Flushed()
	if finalSeq != 3 || cpFinal.StoreSeq != finalSeq || cpFinal.Offset != finalOff {
		t.Fatalf("final cp = %+v, want seq=%d offset=%d", cpFinal, finalSeq, finalOff)
	}
}

func TestCheckpointLocalModeFinalFlush(t *testing.T) {
	s, _, cpFile := newCheckpointTestSession(t, config.ModeLocal, nil)
	t.Cleanup(func() {
		s.stopCheckpoint()
		if s.store != nil {
			_ = s.store.Close()
		}
	})
	s.startCheckpoint()

	commitToIntake(t, s.intake, s.cfg, backfillStart(), backfillScalar(1), backfillScalar(2))
	if err := s.finish(context.Background(), &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("finish: %v", err)
	}

	cp, ok, err := store.LoadCheckpoint(cpFile, s.run.RunID)
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint = (%v, %v), want (nil, true)", err, ok)
	}
	finalSeq, finalOff := s.store.Flushed()
	if finalSeq != 4 || cp.StoreSeq != 4 || cp.Offset != finalOff {
		t.Fatalf("local final cp = %+v, want seq=4 off=%d (flushed=%d)", cp, finalOff, finalSeq)
	}
}

func TestCheckpointWriteFailureWarnsAndRetries(t *testing.T) {
	s, _, cpFile := newCheckpointTestSession(t, config.ModeLocal, nil)
	t.Cleanup(func() {
		_ = s.store.Close()
	})
	// 用同名目录制造 WriteCheckpoint 失败（rename 到非空目录失败）。
	if err := os.MkdirAll(cpFile, 0o755); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	if err := os.WriteFile(filepath.Join(cpFile, "block"), []byte("x"), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	commitToIntake(t, s.intake, s.cfg, backfillStart())
	waitStoreFlushed(t, s.store, 1)

	s.persistCheckpoint()
	if s.checkpointHasWrite || !s.checkpointWarned {
		t.Fatalf("after failed write: hasWrite=%v warned=%v, want (false, true)", s.checkpointHasWrite, s.checkpointWarned)
	}

	// 清除障碍后，下一次 persistCheckpoint 应成功恢复并重置告警标记。
	if err := os.RemoveAll(cpFile); err != nil {
		t.Fatalf("RemoveAll: %v", err)
	}
	s.persistCheckpoint()
	if !s.checkpointHasWrite || s.checkpointWarned {
		t.Fatalf("after recovery: hasWrite=%v warned=%v, want (true, false)", s.checkpointHasWrite, s.checkpointWarned)
	}
	cp, ok, err := store.LoadCheckpoint(cpFile, s.run.RunID)
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint = (%v, %v), want (nil, true)", err, ok)
	}
	if cp.StoreSeq != 1 {
		t.Fatalf("cp.StoreSeq = %d, want 1", cp.StoreSeq)
	}
}

func TestCheckpointSkipStoreProducesNoFile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 1)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-1")

	mgr := newTestManager()
	cfg := startSettings(t, cloud)
	cfg.SkipStore = true

	result, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Project:   "demo",
		Workspace: "alice",
		Id:        "run-1",
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, result.SessionID) })

	s, err := mgr.session(result.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	if s.checkpointFile != "" {
		t.Fatalf("checkpointFile = %q, want empty under skip_store", s.checkpointFile)
	}
	matches, err := filepath.Glob(filepath.Join(cfg.RunDir, "*.checkpoint.json"))
	if err != nil {
		t.Fatalf("Glob: %v", err)
	}
	if len(matches) != 0 {
		t.Fatalf("unexpected checkpoint files under skip_store: %v", matches)
	}
}

func commitToIntake(t *testing.T, in *buffer.Intake, cfg config.Settings, recs ...*recordv1.Record) {
	t.Helper()
	tk, err := in.Reserve(context.Background(), demandFor(recs, cfg))
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := in.Commit(tk, recs); err != nil {
		t.Fatalf("Commit: %v", err)
	}
}
