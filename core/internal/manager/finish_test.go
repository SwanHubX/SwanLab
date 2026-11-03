package manager

import (
	"context"
	"errors"
	"io"
	"os"
	"path/filepath"
	"sync"
	"testing"
	"time"

	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

// collectingSender 记录所有上传批次，模拟一条永不失败的云端链路。
type collectingSender struct {
	mu      sync.Mutex
	records []*recordv1.Record
}

func (c *collectingSender) Upload(_ context.Context, records []*recordv1.Record) uploader.Result {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.records = append(c.records, records...)
	return uploader.Result{Confirmed: records}
}

func (c *collectingSender) snapshot() []*recordv1.Record {
	c.mu.Lock()
	defer c.mu.Unlock()
	return append([]*recordv1.Record(nil), c.records...)
}

type finishOptions struct {
	skipStore   bool
	uploadQueue int
}

type finishHarness struct {
	mgr     *Manager
	id      string
	session *Session
	runFile string
	writer  *store.Writer
	uploads *collectingSender
}

// newFinishHarness 装配真实 store/intake/uploader，用于验证收尾序列的顺序与刷盘语义。
func newFinishHarness(t *testing.T, opts finishOptions) *finishHarness {
	t.Helper()
	if opts.uploadQueue <= 0 {
		opts.uploadQueue = 16
	}
	cfg := config.Defaults()
	cfg.RunID = "finish"
	cfg.RunDir = t.TempDir()
	cfg.RecordBatch = 8
	cfg.RecordInterval = 0.01
	cfg.SkipStore = opts.skipStore

	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	runFile := ""

	var writer *store.Writer
	var async *store.AsyncWriter
	if !cfg.SkipStore {
		runFile, err = run.RunFile()
		if err != nil {
			t.Fatalf("RunFile: %v", err)
		}
		if err = os.MkdirAll(run.RunDir, 0o755); err != nil {
			t.Fatalf("MkdirAll: %v", err)
		}
		writer, err = store.Create(runFile)
		if err != nil {
			t.Fatalf("store.Create: %v", err)
		}
		async, err = store.NewAsyncWriter(writer, store.AsyncOptions{Queue: 16})
		if err != nil {
			t.Fatalf("store.NewAsyncWriter: %v", err)
		}
	}

	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.run = run
	s.metrics = metrics.New()
	s.counter = counter.New(0)
	s.epoch = counter.New(0)
	s.watcher = watcher.NullWatcher{}
	s.tracker = uploader.NewTracker()

	b := budget.New(budget.Limits{})
	var storeIn chan<- store.Request
	var flushed func() uint64
	if async != nil {
		storeIn = async.In()
		flushed = func() uint64 {
			seq, _ := async.Flushed()
			return seq
		}
	}
	intake, err := buffer.New(buffer.Options{
		SkipStore:   cfg.SkipStore,
		Store:       storeIn,
		UploadQueue: opts.uploadQueue,
		Budget:      b,
		Flushed:     flushed,
		OnAdmit:     s.tracker.AddTotalRecords,
		QueueWait:   20 * time.Millisecond,
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}

	uploads := &collectingSender{}
	senders := map[record.Type]uploader.Sender{
		record.TypeColumn: uploads,
		record.TypeScalar: uploads,
		record.TypeLog:    uploads,
		record.TypeMedia:  uploads,
		record.TypeSave:   uploads,
	}
	up, err := uploader.New(uploader.Options{
		In:       intake.UploadOut(),
		Budget:   b,
		Senders:  senders,
		Batch:    4,
		Interval: 10 * time.Millisecond,
		Tracker:  s.tracker,
		Context:  s.ctx,
		OnDrained: func() {
			s.markDrained()
		},
	})
	if err != nil {
		t.Fatalf("uploader.New: %v", err)
	}

	s.store = async
	s.intake = intake
	s.uploader = up
	s.drainFn = up.Wait
	s.runFile = runFile
	if async != nil {
		s.startRecovery()
	}

	mgr := newTestManager()
	id, err := mgr.install(s)
	if err != nil {
		t.Fatalf("install: %v", err)
	}
	t.Cleanup(func() {
		s.stopRecovery()
		s.cancel()
		_ = up.Close()
		if async != nil {
			_ = async.Close()
		}
	})
	return &finishHarness{mgr: mgr, id: id, session: s, runFile: runFile, writer: writer, uploads: uploads}
}

func waitDrained(t *testing.T, s *Session) {
	t.Helper()
	select {
	case <-s.drained:
	case <-time.After(3 * time.Second):
		t.Fatal("session did not drain")
	}
}

func readStore(t *testing.T, path string) []*recordv1.Record {
	t.Helper()
	reader, err := store.Open(path)
	if err != nil {
		t.Fatalf("store.Open: %v", err)
	}
	defer func() { _ = reader.Close() }()
	var records []*recordv1.Record
	for {
		payload, err := reader.Scan()
		if errors.Is(err, io.EOF) {
			return records
		}
		if err != nil {
			t.Fatalf("scan: %v", err)
		}
		rec := &recordv1.Record{}
		if err := proto.Unmarshal(payload, rec); err != nil {
			t.Fatalf("unmarshal record: %v", err)
		}
		records = append(records, rec)
	}
}

func finishRecord(state runv1.RunState, message string) *runv1.FinishRecord {
	return &runv1.FinishRecord{
		State:      state,
		Error:      message,
		FinishedAt: timestamppb.New(time.Unix(1_700_000_000, 0)),
	}
}

func TestDeliverRunFinishFlushesAndIsIdempotent(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("finish/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}

	rec := finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")
	if err := h.mgr.DeliverRunFinish(ctx, h.id, rec); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	// 重复投递必须幂等，且不能产生第二条 FinishRecord。
	if err := h.mgr.DeliverRunFinish(ctx, h.id, rec); err != nil {
		t.Fatalf("repeated DeliverRunFinish: %v", err)
	}
	waitDrained(t, h.session)

	records := readStore(t, h.runFile)
	if len(records) == 0 {
		t.Fatal("run file is empty")
	}
	last := records[len(records)-1]
	if last.GetNum() != record.NumFinish {
		t.Fatalf("last record num = %d, want %d", last.GetNum(), record.NumFinish)
	}
	if last.GetFinish().GetState() != runv1.RunState_RUN_STATE_FINISHED {
		t.Fatalf("finish state = %v, want FINISHED", last.GetFinish().GetState())
	}
	sawFinish := 0
	for _, item := range records {
		if item.GetNum() == record.NumFinish {
			sawFinish++
		}
	}
	if sawFinish != 1 {
		t.Fatalf("stored %d finish records, want 1", sawFinish)
	}

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("finish/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); !errors.Is(err, ErrNotRunning) {
		t.Fatalf("upsert after finish err=%v, want ErrNotRunning", err)
	}
}

func TestDeliverRunFinishUploadsErrorLogBeforeDrain(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	rec := finishRecord(runv1.RunState_RUN_STATE_CRASHED, "boom")
	if err := h.mgr.DeliverRunFinish(ctx, h.id, rec); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	waitDrained(t, h.session)

	uploads := h.uploads.snapshot()
	if len(uploads) != 1 || uploads[0].GetLog() == nil {
		t.Fatalf("uploads = %v, want one error log", uploads)
	}
	log := uploads[0].GetLog()
	if log.GetLine() != "boom" {
		t.Fatalf("error log line = %q, want boom", log.GetLine())
	}
	if log.GetLevel() != terminalv1.LogLevel_LOG_LEVEL_ERROR {
		t.Fatalf("error log level = %v, want ERROR", log.GetLevel())
	}
	if got := log.GetTimestamp().AsTime(); !got.Equal(rec.GetFinishedAt().AsTime()) {
		t.Fatalf("error log timestamp = %v, want %v", got, rec.GetFinishedAt().AsTime())
	}

	records := readStore(t, h.runFile)
	if len(records) != 1 || records[0].GetFinish().GetState() != runv1.RunState_RUN_STATE_CRASHED {
		t.Fatalf("stored records = %v, want one CRASHED finish record", records)
	}
}

func TestDeliverRunFinishUsesDefaultErrorText(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_ABORTED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	waitDrained(t, h.session)

	uploads := h.uploads.snapshot()
	if len(uploads) != 1 || uploads[0].GetLog().GetLine() != finishUnknownError {
		t.Fatalf("error log = %v, want default message %q", uploads, finishUnknownError)
	}
}

func TestDeliverRunFinishStoreFailureDegradesOnly(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	// 模拟磁盘持续写入失败：AsyncWriter 只能降级，收尾仍须完成上传侧排空。
	if err := h.writer.Close(); err != nil {
		t.Fatalf("close writer: %v", err)
	}
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish with broken store: %v", err)
	}
	waitDrained(t, h.session)
	if !h.session.storeDegraded.Load() {
		t.Fatal("store failure must set storeDegraded")
	}
	if state := h.session.stateSnapshot(); state != RunStateDrained {
		t.Fatalf("state = %s, want drained", state)
	}
}

func TestDeliverRunFinishIntakeCloseFailureMarksFailed(t *testing.T) {
	h := newFinishHarness(t, finishOptions{uploadQueue: 2})
	// 先停掉 uploader，让上传通道不再被消费；随后填满唯一槽位。
	if err := h.session.uploader.Close(); err != nil && !errors.Is(err, uploader.ErrNotDrained) {
		t.Fatalf("close uploader: %v", err)
	}
	setupCtx := context.Background()
	if err := h.mgr.UpsertScalars(setupCtx, h.id, []*datav1.ScalarRecord{
		scalarRecord("finish/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	endPath := filepath.Join(t.TempDir(), "end.txt")
	if err := os.WriteFile(endPath, []byte("end"), 0o644); err != nil {
		t.Fatalf("write end file: %v", err)
	}
	if err := h.mgr.UpsertSaves(setupCtx, h.id, []*savev1.SaveRecord{{
		Name:       "end.txt",
		SourcePath: endPath,
		Policy:     savev1.SavePolicy_SAVE_POLICY_END,
		Type:       savev1.SaveType_SAVE_TYPE_CUSTOM,
	}}); err != nil {
		t.Fatalf("UpsertSaves: %v", err)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 50*time.Millisecond)
	defer cancel()
	err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
	if err == nil || !errors.Is(err, ErrRunFailed) {
		t.Fatalf("DeliverRunFinish err=%v, want ErrRunFailed", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateFailed {
		t.Fatalf("state = %s, want failed", state)
	}
}

func TestDeliverRunFinishConcurrentCallsWriteOnce(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	const callers = 4
	var wg sync.WaitGroup
	errs := make(chan error, callers)
	for range callers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			errs <- h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
		}()
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		if err != nil {
			t.Fatalf("concurrent DeliverRunFinish: %v", err)
		}
	}
	waitDrained(t, h.session)

	finishCount := 0
	for _, item := range readStore(t, h.runFile) {
		if item.GetNum() == record.NumFinish {
			finishCount++
		}
	}
	if finishCount != 1 {
		t.Fatalf("stored %d finish records, want 1", finishCount)
	}
}

func TestDeliverRunFinishRejectsNilRecord(t *testing.T) {
	m := newTestManager()
	handle := installBare(t, m, sessionOptions{})
	if err := m.DeliverRunFinish(context.Background(), handle, nil); !errors.Is(err, ErrInvalidFinishRecord) {
		t.Fatalf("nil finish record err=%v, want ErrInvalidFinishRecord", err)
	}
	if err := m.DeliverRunFinish(context.Background(), "", finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); !errors.Is(err, ErrEmptySessionID) {
		t.Fatalf("empty session err=%v, want ErrEmptySessionID", err)
	}
}

func TestDeliverRunFinishSkipStore(t *testing.T) {
	h := newFinishHarness(t, finishOptions{skipStore: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertSaves(ctx, h.id, []*savev1.SaveRecord{{
		Name:       "config.yaml",
		Policy:     savev1.SavePolicy_SAVE_POLICY_NOW,
		Type:       savev1.SaveType_SAVE_TYPE_CONFIG,
		Payload:    []byte("key: value\n"),
		SourcePath: "",
	}}); err != nil {
		t.Fatalf("UpsertSaves: %v", err)
	}
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	waitDrained(t, h.session)
	if state := h.session.stateSnapshot(); state != RunStateDrained {
		t.Fatalf("state = %s, want drained", state)
	}
	if got := len(h.uploads.snapshot()); got != 1 {
		t.Fatalf("uploaded %d records, want only the config save", got)
	}
}

// Shutdown 面向进程退出：排空会话、刷盘且不写 FinishRecord，也不上报实验终态。
func TestManagerShutdownDrainsActiveSessions(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("shutdown/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	var stops int
	h.session.releaseFn = func(context.Context) error {
		stops++
		return nil
	}
	if err := h.mgr.Shutdown(ctx); err != nil {
		t.Fatalf("Shutdown: %v", err)
	}
	if stops != 0 {
		t.Fatalf("Shutdown must not call stop_experiment, got %d call(s)", stops)
	}
	if state := h.session.stateSnapshot(); state != RunStateDrained {
		t.Fatalf("state = %s, want drained", state)
	}

	records := readStore(t, h.runFile)
	if len(records) == 0 {
		t.Fatal("shutdown must flush accepted records to the local store")
	}
	for _, item := range records {
		if item.GetNum() == record.NumFinish {
			t.Fatal("shutdown must not synthesize a FinishRecord")
		}
	}
	if got := len(h.uploads.snapshot()); got == 0 {
		t.Fatal("shutdown must drain accepted records through the uploader")
	}
}
