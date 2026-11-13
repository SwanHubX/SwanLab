package manager

import (
	"context"
	"errors"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"sync"
	"testing"
	"time"

	"go.uber.org/goleak"
	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

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
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

func newTestManager() *Manager { return New(Options{}) }

func installBare(t *testing.T, m *Manager, opts sessionOptions) string {
	t.Helper()
	handle, err := m.install(newSession(opts))
	if err != nil {
		t.Fatalf("install session: %v", err)
	}
	return handle
}

func serveAuth(cloud *clienttest.Cloud) {
	cloud.ServeJSON(http.MethodGet, "/api/auth/verify", http.StatusOK, map[string]any{
		"uid":       1,
		"username":  "alice",
		"createdAt": "2026-01-01T00:00:00Z",
	})
}

func startSettings(t *testing.T, cloud *clienttest.Cloud) config.Settings {
	t.Helper()
	settings := config.Defaults()
	settings.RunID = "run-1"
	settings.RunDir = filepath.Join(t.TempDir(), "run")
	settings.RecordBatch = 16
	settings.RecordInterval = 0.01
	settings.APIKey = "smoke-key"
	settings.APIHost = cloud.URL()
	settings.WebHost = cloud.URL()
	return settings
}

func runFilePath(t *testing.T, settings config.Settings) string {
	t.Helper()
	run, err := runctx.New(settings)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	path, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}
	return path
}

// confirmAll 立即确认整批记录，用于构造空闲 uploader。
type confirmAll struct{}

func (confirmAll) Upload(_ context.Context, records []*recordv1.Record) uploader.Result {
	return uploader.Result{Confirmed: records}
}

func newIdleUploader(t *testing.T) *uploader.Uploader {
	t.Helper()
	up, err := uploader.New(uploader.Options{
		In:       make(chan *recordv1.Record),
		Budget:   budget.New(budget.Limits{}),
		Senders:  map[record.Type]uploader.Sender{record.TypeScalar: confirmAll{}},
		Interval: time.Hour,
	})
	if err != nil {
		t.Fatalf("uploader.New: %v", err)
	}
	t.Cleanup(func() { _ = up.Close() })
	return up
}

func readStore(t *testing.T, path string) []*recordv1.Record {
	t.Helper()
	reader, err := store.Open(path)
	if err != nil {
		t.Fatalf("store.Open(%q): %v", path, err)
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

func scalarRecord(key string, step int64, typ columnv1.ColumnType) *datav1.ScalarRecord {
	return &datav1.ScalarRecord{Key: key, Step: step, Type: typ}
}

func finishRecord(state runv1.RunState, message string) *runv1.FinishRecord {
	return &runv1.FinishRecord{
		State:      state,
		Error:      message,
		FinishedAt: timestamppb.New(time.Unix(1_700_000_000, 0)),
	}
}

// collectingSender 记录所有上传批次，模拟一条永不失败的云端链路。
type collectingSender struct {
	mu      sync.Mutex
	records []*recordv1.Record
	notify  chan struct{}
}

func (c *collectingSender) Upload(_ context.Context, records []*recordv1.Record) uploader.Result {
	c.mu.Lock()
	c.records = append(c.records, records...)
	c.mu.Unlock()
	select {
	case c.notify <- struct{}{}:
	default:
	}
	return uploader.Result{Confirmed: records}
}

// waitFor 等到至少 n 条记录送达；用发送信号同步，不做轮询。
func (c *collectingSender) waitFor(t *testing.T, n int) {
	t.Helper()
	deadline := time.NewTimer(5 * time.Second)
	defer deadline.Stop()
	for len(c.snapshot()) < n {
		select {
		case <-c.notify:
		case <-deadline.C:
			t.Fatalf("uploader delivered %d records, want at least %d", len(c.snapshot()), n)
		}
	}
}

func (c *collectingSender) snapshot() []*recordv1.Record {
	c.mu.Lock()
	defer c.mu.Unlock()
	return append([]*recordv1.Record(nil), c.records...)
}

// gatedSender 在首次上传时阻塞，用于制造"上传通道持续满"的瞬时故障。
type gatedSender struct {
	inner   uploader.Sender
	gate    chan struct{}
	entered chan struct{}
	once    sync.Once
}

func (g *gatedSender) Upload(ctx context.Context, records []*recordv1.Record) uploader.Result {
	g.once.Do(func() { close(g.entered) })
	select {
	case <-g.gate:
	case <-ctx.Done():
		return uploader.Result{}
	}
	return g.inner.Upload(ctx, records)
}

type finishOptions struct {
	skipStore     bool
	mode          config.Mode
	uploadQueue   int
	queueWait     time.Duration
	blockUpload   bool
	budgetRecords int64
}

type finishHarness struct {
	mgr           *Manager
	id            string
	session       *Session
	runFile       string
	writer        *store.Writer
	uploads       *collectingSender
	gate          chan struct{}
	uploadEntered chan struct{}
}

// newFinishHarness 装配真实 store/intake/uploader，用于验证收尾序列的顺序与刷盘语义。
func newFinishHarness(t *testing.T, opts finishOptions) *finishHarness {
	t.Helper()
	if opts.uploadQueue <= 0 {
		opts.uploadQueue = 16
	}
	if opts.queueWait <= 0 {
		opts.queueWait = 20 * time.Millisecond
	}
	cfg := config.Defaults()
	cfg.RunID = "finish"
	cfg.RunDir = t.TempDir()
	cfg.RecordBatch = 8
	cfg.RecordInterval = 0.01
	cfg.SkipStore = opts.skipStore
	cfg.Mode = opts.mode
	storeOnly := cfg.NormalizedMode() != config.ModeOnline

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

	b := budget.New(budget.Limits{Records: opts.budgetRecords})
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
		StoreOnly:   storeOnly,
		Store:       storeIn,
		UploadQueue: opts.uploadQueue,
		Budget:      b,
		Flushed:     flushed,
		OnAdmit:     s.tracker.AddTotalRecords,
		QueueWait:   opts.queueWait,
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}

	uploads := &collectingSender{notify: make(chan struct{}, 1)}
	var gate chan struct{}
	var entered chan struct{}
	var up *uploader.Uploader
	if !storeOnly {
		var sender uploader.Sender = uploads
		if opts.blockUpload {
			gate = make(chan struct{})
			entered = make(chan struct{})
			sender = &gatedSender{inner: uploads, gate: gate, entered: entered}
		}
		senders := map[record.Type]uploader.Sender{
			record.TypeColumn: sender,
			record.TypeScalar: sender,
			record.TypeLog:    sender,
			record.TypeMedia:  sender,
			record.TypeSave:   sender,
		}
		up, err = uploader.New(uploader.Options{
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
	}

	s.store = async
	s.intake = intake
	s.uploader = up
	if up != nil {
		s.drainFn = up.Wait
	} else {
		// local/offline 没有 uploader，收尾完成即排空。
		s.drainFn = func(context.Context) error {
			s.markDrained()
			return nil
		}
	}
	s.runFile = runFile
	if async != nil && !storeOnly {
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
		if up != nil {
			_ = up.Close()
		}
		if async != nil {
			_ = async.Close()
		}
	})
	return &finishHarness{
		mgr: mgr, id: id, session: s, runFile: runFile, writer: writer,
		uploads: uploads, gate: gate, uploadEntered: entered,
	}
}
