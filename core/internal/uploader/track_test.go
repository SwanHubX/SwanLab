package uploader

import (
	"context"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// senderFunc 是 Sender 的函数适配器。
type senderFunc func(context.Context, []*recordv1.Record) Result

func (f senderFunc) Upload(ctx context.Context, records []*recordv1.Record) Result {
	return f(ctx, records)
}

// scriptedSender 按调用序号返回预置结果，超出脚本后重复最后一条。
type scriptedSender struct {
	mu      sync.Mutex
	scripts []func([]*recordv1.Record) Result
	calls   [][]*recordv1.Record
}

func (s *scriptedSender) Upload(_ context.Context, records []*recordv1.Record) Result {
	s.mu.Lock()
	index := len(s.calls)
	s.calls = append(s.calls, records)
	var step func([]*recordv1.Record) Result
	if index < len(s.scripts) {
		step = s.scripts[index]
	} else if len(s.scripts) > 0 {
		step = s.scripts[len(s.scripts)-1]
	}
	s.mu.Unlock()

	if step == nil {
		return Result{Confirmed: records}
	}
	return step(records)
}

// quietThrottle 返回不输出任何内容的告警节流器。
func quietThrottle() *WarningThrottle {
	throttle := newWarningThrottle(time.Now, time.Hour)
	throttle.warnf = func(string) {}
	throttle.infof = func(string) {}
	return throttle
}

// retryNoJitter 返回无抖动的退避：每次失败固定等待 base。
func retryNoJitter(base time.Duration) RetryPolicy {
	return RetryPolicy{Base: base, Multiplier: 2, Max: base, JitterLow: 1, JitterHigh: 1}
}

type regFixture struct {
	reg      *Registry
	ledger   *budget.Budget
	tracker  *Tracker
	progress chan struct{}
}

func newRegFixture(t *testing.T, cfg registryConfig) *regFixture {
	t.Helper()
	if cfg.budget == nil {
		cfg.budget = budget.New(budget.Limits{})
	}
	if cfg.tracker == nil {
		cfg.tracker = NewTracker()
	}
	if cfg.warn == nil {
		cfg.warn = quietThrottle()
	}
	if cfg.ctx == nil {
		cfg.ctx = context.Background()
	}
	progress := make(chan struct{}, 1)
	cfg.notify = func() {
		select {
		case progress <- struct{}{}:
		default:
		}
	}
	reg := newRegistry(cfg)
	t.Cleanup(func() {
		reg.stop()
		reg.wait()
	})
	return &regFixture{reg: reg, ledger: cfg.budget, tracker: cfg.tracker, progress: progress}
}

// own 模拟 intake 为记录预留配额。
func (f *regFixture) own(t *testing.T, records ...*recordv1.Record) {
	t.Helper()
	for _, rec := range records {
		if err := f.ledger.Acquire(context.Background(), 1, record.Size(rec)); err != nil {
			t.Fatalf("acquire: %v", err)
		}
	}
}

// takeOne 取走某条轨道当前可出队的批次。
func (f *regFixture) takeOne(t *testing.T, typ record.Type) []*recordv1.Record {
	t.Helper()
	track := f.reg.tracks[typ]
	batch, _, exit := f.reg.take(track)
	if exit || len(batch) == 0 {
		t.Fatalf("expected a batch from %s, got exit=%v len=%d", typ, exit, len(batch))
	}
	return batch
}

// waitDrained 等待 Registry 判定排空。
func (f *regFixture) waitDrained(t *testing.T) {
	t.Helper()
	timeout := time.After(2 * time.Second)
	for {
		if f.reg.drained() {
			return
		}
		select {
		case <-f.progress:
		case <-timeout:
			t.Fatal("timeout waiting for drained")
		}
	}
}

func waitSignal(t *testing.T, ch <-chan struct{}, what string) {
	t.Helper()
	select {
	case <-ch:
	case <-time.After(2 * time.Second):
		t.Fatalf("timeout waiting for %s", what)
	}
}

func TestRegistryDeterministicTrackOrder(t *testing.T) {
	senders := map[record.Type]Sender{}
	for _, typ := range dispatchOrder {
		senders[typ] = senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
			return Result{Confirmed: records}
		})
	}
	f := newRegFixture(t, registryConfig{senders: senders, batch: 10})

	stats := f.reg.stats()
	if len(stats) != len(dispatchOrder) {
		t.Fatalf("tracks = %d, want %d", len(stats), len(dispatchOrder))
	}
	for i, typ := range dispatchOrder {
		if stats[i].Type != typ {
			t.Fatalf("track[%d] = %s, want %s", i, stats[i].Type, typ)
		}
	}
}

func TestRegistryUnregisteredSenderIgnored(t *testing.T) {
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeStart: &scriptedSender{}},
		batch:   10,
	})
	if got := len(f.reg.stats()); got != 0 {
		t.Fatalf("tracks = %d, want 0", got)
	}
}

func TestRegistryDrainedRequiresIngestClosed(t *testing.T) {
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: &scriptedSender{}},
		batch:   10,
	})

	if f.reg.drained() {
		t.Fatal("must not report drained before ingest closes")
	}
	f.reg.markIngestClosed()
	if !f.reg.drained() {
		t.Fatal("empty registry with closed ingest should be drained")
	}

	recs := []*recordv1.Record{newScalar(1, "a")}
	f.own(t, recs...)
	f.reg.dispatch([]recordGroup{{typ: record.TypeScalar, records: recs}})
	if f.reg.drained() {
		t.Fatal("queued records must keep registry undrained")
	}
}

func TestUploadAccountsPartialSuccess(t *testing.T) {
	clock := newFakeClock()
	sender := &scriptedSender{scripts: []func([]*recordv1.Record) Result{
		func(records []*recordv1.Record) Result {
			return Result{Confirmed: []*recordv1.Record{records[0], records[2]}, Remaining: []*recordv1.Record{records[1]}}
		},
		func(records []*recordv1.Record) Result { return Result{Confirmed: records} },
	}}
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: sender},
		batch:   10,
		retry:   retryNoJitter(time.Second),
		now:     clock.now,
	})

	recs := []*recordv1.Record{newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c")}
	f.own(t, recs...)
	f.reg.dispatch([]recordGroup{{typ: record.TypeScalar, records: recs}})
	batch := f.takeOne(t, record.TypeScalar)

	f.reg.upload(f.reg.tracks[record.TypeScalar], batch)

	if used, _ := f.ledger.Usage(); used != 1 {
		t.Fatalf("budget usage = %d, want 1 (only the retried record)", used)
	}
	if got := f.tracker.Snapshot().UploadedRecords; got != 2 {
		t.Fatalf("uploaded records = %d, want 2", got)
	}
	stats := f.reg.stats()[0]
	if stats.Queue != 1 || stats.Confirmed != 2 || stats.FailedRounds != 1 {
		t.Fatalf("stats = %+v, want queue=1 confirmed=2 failed=1", stats)
	}
	if want := clock.now().Add(time.Second); !stats.RetryAt.Equal(want) {
		t.Fatalf("retry at = %v, want %v", stats.RetryAt, want)
	}
	if pending := f.reg.pendingRecords(); len(pending) != 1 || pending[0] != recs[1] {
		t.Fatalf("pending = %v, want the remaining record", pending)
	}

	// 退避窗口未到不能取批；推进时钟后重试成功，配额归零、失败轮次复位。
	if _, _, exit := f.reg.take(f.reg.tracks[record.TypeScalar]); exit {
		t.Fatal("unexpected worker exit")
	}
	if got := f.reg.tracks[record.TypeScalar].inflight; got != 0 {
		t.Fatalf("inflight = %d, want 0 while backing off", got)
	}
	clock.advance(time.Second)
	batch = f.takeOne(t, record.TypeScalar)
	f.reg.upload(f.reg.tracks[record.TypeScalar], batch)

	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0 after retry", used)
	}
	stats = f.reg.stats()[0]
	if stats.Queue != 0 || stats.Confirmed != 3 || stats.FailedRounds != 0 || !stats.RetryAt.IsZero() {
		t.Fatalf("stats after retry = %+v", stats)
	}
}

func TestUploadRejectsInconsistentResult(t *testing.T) {
	sender := &scriptedSender{scripts: []func([]*recordv1.Record) Result{
		func(records []*recordv1.Record) Result {
			return Result{Confirmed: records[:2]} // 少算一条，守恒校验必须拦截
		},
		func(records []*recordv1.Record) Result { return Result{Confirmed: records} },
	}}
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: sender},
		batch:   10,
		retry:   retryNoJitter(time.Millisecond),
	})

	recs := []*recordv1.Record{newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c")}
	f.own(t, recs...)
	f.reg.dispatch([]recordGroup{{typ: record.TypeScalar, records: recs}})
	f.reg.upload(f.reg.tracks[record.TypeScalar], f.takeOne(t, record.TypeScalar))

	if used, _ := f.ledger.Usage(); used != 3 {
		t.Fatalf("budget usage = %d, want 3 (nothing confirmed)", used)
	}
	if got := f.reg.stats()[0].Queue; got != 3 {
		t.Fatalf("queue = %d, want whole batch back", got)
	}
	failures := f.reg.lastFailures()
	if len(failures) != 1 || failures[0].Message != "sender result does not account for every record" {
		t.Fatalf("failures = %+v, want conservation failure", failures)
	}

	// 下一轮正常确认后全部归还。
	f.reg.tracks[record.TypeScalar].retryAt = time.Time{}
	f.reg.upload(f.reg.tracks[record.TypeScalar], f.takeOne(t, record.TypeScalar))
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0", used)
	}
}

func TestUploadRecoversSenderPanic(t *testing.T) {
	var attempts atomic.Int32
	sender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		if attempts.Add(1) == 1 {
			panic("boom")
		}
		return Result{Confirmed: records}
	})
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: sender},
		batch:   10,
		retry:   retryNoJitter(time.Millisecond),
	})

	recs := []*recordv1.Record{newScalar(1, "a"), newScalar(2, "b")}
	f.own(t, recs...)
	f.reg.dispatch([]recordGroup{{typ: record.TypeScalar, records: recs}})
	f.reg.upload(f.reg.tracks[record.TypeScalar], f.takeOne(t, record.TypeScalar))

	if used, _ := f.ledger.Usage(); used != 2 {
		t.Fatalf("budget usage = %d, want 2 after panic", used)
	}
	failures := f.reg.lastFailures()
	if len(failures) != 1 || failures[0].Message != "sender panic: boom" {
		t.Fatalf("failures = %+v, want panic failure", failures)
	}
	if got := f.reg.stats()[0].Queue; got != 2 {
		t.Fatalf("queue = %d, want whole batch back", got)
	}

	f.reg.tracks[record.TypeScalar].retryAt = time.Time{}
	f.reg.upload(f.reg.tracks[record.TypeScalar], f.takeOne(t, record.TypeScalar))
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0 after recovery", used)
	}
}

func TestRemainingKeepsBatchOrder(t *testing.T) {
	sender := &scriptedSender{scripts: []func([]*recordv1.Record) Result{
		func(records []*recordv1.Record) Result {
			// 故意倒序返回 remaining
			return Result{Confirmed: []*recordv1.Record{records[1]}, Remaining: []*recordv1.Record{records[2], records[0]}}
		},
		func(records []*recordv1.Record) Result { return Result{Confirmed: records} },
	}}
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: sender},
		batch:   10,
		retry:   retryNoJitter(time.Millisecond),
	})

	recs := []*recordv1.Record{newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c")}
	f.own(t, recs...)
	f.reg.dispatch([]recordGroup{{typ: record.TypeScalar, records: recs}})
	f.reg.upload(f.reg.tracks[record.TypeScalar], f.takeOne(t, record.TypeScalar))

	pending := f.reg.pendingRecords()
	if len(pending) != 2 || pending[0] != recs[0] || pending[1] != recs[2] {
		t.Fatalf("pending order = %v, want [a, c]", pending)
	}
}

func TestSuspendResumeFreezesRetryAt(t *testing.T) {
	clock := newFakeClock()
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: &scriptedSender{}},
		batch:   10,
		now:     clock.now,
	})
	track := f.reg.tracks[record.TypeScalar]
	start := clock.now()
	track.retryAt = start.Add(100 * time.Millisecond)

	f.reg.suspend()
	clock.advance(40 * time.Millisecond)
	f.reg.resume()

	if want := start.Add(140 * time.Millisecond); !track.retryAt.Equal(want) {
		t.Fatalf("retry at = %v, want %v (frozen during suspend)", track.retryAt, want)
	}

	// 重复 resume 不产生第二次顺延。
	f.reg.resume()
	if want := start.Add(140 * time.Millisecond); !track.retryAt.Equal(want) {
		t.Fatalf("retry at after duplicate resume = %v, want %v", track.retryAt, want)
	}
}

func TestRegistryFIFOBatches(t *testing.T) {
	batches := make(chan []int64, 4)
	sender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		nums := make([]int64, 0, len(records))
		for _, rec := range records {
			nums = append(nums, rec.GetNum())
		}
		batches <- nums
		return Result{Confirmed: records}
	})
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: sender},
		batch:   2,
	})

	recs := []*recordv1.Record{
		newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c"), newScalar(4, "d"), newScalar(5, "e"),
	}
	f.own(t, recs...)
	f.reg.start()
	f.reg.dispatch([]recordGroup{{typ: record.TypeScalar, records: recs}})
	f.reg.markIngestClosed()

	want := [][]int64{{1, 2}, {3, 4}, {5}}
	for _, expected := range want {
		select {
		case got := <-batches:
			if len(got) != len(expected) {
				t.Fatalf("batch = %v, want %v", got, expected)
			}
			for i := range expected {
				if got[i] != expected[i] {
					t.Fatalf("batch = %v, want %v", got, expected)
				}
			}
		case <-time.After(2 * time.Second):
			t.Fatalf("timeout waiting for batch %v", expected)
		}
	}

	f.waitDrained(t)
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0", used)
	}
}

func TestColumnGateBlocksDataTracks(t *testing.T) {
	release := make(chan struct{})
	columnStarted := make(chan struct{}, 1)
	order := make(chan record.Type, 4)

	column := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		columnStarted <- struct{}{}
		<-release
		order <- record.TypeColumn
		return Result{Confirmed: records}
	})
	scalar := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		order <- record.TypeScalar
		return Result{Confirmed: records}
	})
	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{
			record.TypeColumn: column,
			record.TypeScalar: scalar,
		},
		batch: 10,
	})

	columnRecs := []*recordv1.Record{newColumn(1, "loss")}
	scalarRecs := []*recordv1.Record{newScalar(2, "loss")}
	f.own(t, columnRecs[0], scalarRecs[0])
	f.reg.start()
	f.reg.dispatch([]recordGroup{
		{typ: record.TypeColumn, records: columnRecs},
		{typ: record.TypeScalar, records: scalarRecs},
	})

	waitSignal(t, columnStarted, "column upload")
	select {
	case typ := <-order:
		t.Fatalf("%s track ran while column was in flight", typ)
	default:
	}

	close(release)
	if got := <-order; got != record.TypeColumn {
		t.Fatalf("first completed track = %s, want column", got)
	}
	if got := <-order; got != record.TypeScalar {
		t.Fatalf("second completed track = %s, want scalar", got)
	}

	f.reg.markIngestClosed()
	f.waitDrained(t)
}

func TestTrackFailureDoesNotBlockOthers(t *testing.T) {
	scalarAttempts := make(chan struct{}, 4)
	var scalarCalls atomic.Int32
	scalar := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		if scalarCalls.Add(1) == 1 {
			scalarAttempts <- struct{}{}
			return Result{Remaining: records}
		}
		scalarAttempts <- struct{}{}
		return Result{Confirmed: records}
	})
	logDone := make(chan struct{}, 1)
	logSender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		logDone <- struct{}{}
		return Result{Confirmed: records}
	})

	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{
			record.TypeScalar: scalar,
			record.TypeLog:    logSender,
		},
		batch: 10,
		retry: retryNoJitter(time.Hour), // 退避期间只能靠 beginFinish 打断
	})

	scalarRecs := []*recordv1.Record{newScalar(1, "loss")}
	logRecs := []*recordv1.Record{newLog(2, "line")}
	f.own(t, scalarRecs[0], logRecs[0])
	f.reg.start()
	f.reg.dispatch([]recordGroup{
		{typ: record.TypeScalar, records: scalarRecs},
		{typ: record.TypeLog, records: logRecs},
	})

	waitSignal(t, scalarAttempts, "first scalar attempt")
	waitSignal(t, logDone, "log upload while scalar backs off")
	if got := scalarCalls.Load(); got != 1 {
		t.Fatalf("scalar attempts = %d, want 1 while log proceeds", got)
	}

	f.reg.beginFinish()
	waitSignal(t, scalarAttempts, "scalar retry after finish")

	f.reg.markIngestClosed()
	f.waitDrained(t)
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0", used)
	}
}

func TestNoRecordsLostAcrossRetries(t *testing.T) {
	const total = 40
	recs := make([]*recordv1.Record, 0, total)
	for i := 0; i < total; i++ {
		recs = append(recs, newScalar(int64(i+1), "k"))
	}

	var mu sync.Mutex
	confirmed := make(map[*recordv1.Record]int, total)
	attempt := 0
	sender := senderFunc(func(_ context.Context, batch []*recordv1.Record) Result {
		mu.Lock()
		attempt++
		first := attempt == 1
		mu.Unlock()

		var done, retry []*recordv1.Record
		for i, rec := range batch {
			if first && i%2 == 1 {
				retry = append(retry, rec)
				continue
			}
			done = append(done, rec)
		}
		mu.Lock()
		for _, rec := range done {
			confirmed[rec]++
		}
		mu.Unlock()
		return Result{Confirmed: done, Remaining: retry}
	})

	f := newRegFixture(t, registryConfig{
		senders: map[record.Type]Sender{record.TypeScalar: sender},
		batch:   7,
		retry:   retryNoJitter(time.Millisecond),
	})
	f.own(t, recs...)
	f.reg.start()
	f.reg.dispatch([]recordGroup{{typ: record.TypeScalar, records: recs}})
	f.reg.markIngestClosed()
	f.waitDrained(t)

	mu.Lock()
	defer mu.Unlock()
	if len(confirmed) != total {
		t.Fatalf("confirmed %d distinct records, want %d", len(confirmed), total)
	}
	for i, rec := range recs {
		if confirmed[rec] != 1 {
			t.Fatalf("record %d confirmed %d times, want exactly once", i, confirmed[rec])
		}
	}
	if got := f.tracker.Snapshot().UploadedRecords; got != total {
		t.Fatalf("uploaded records = %d, want %d", got, total)
	}
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0", used)
	}
}

func TestConserved(t *testing.T) {
	a, b := newScalar(1, "a"), newScalar(2, "b")
	batch := []*recordv1.Record{a, b}

	cases := []struct {
		name   string
		result Result
		want   bool
	}{
		{"all confirmed", Result{Confirmed: batch}, true},
		{"all remaining", Result{Remaining: batch}, true},
		{"split", Result{Confirmed: []*recordv1.Record{a}, Remaining: []*recordv1.Record{b}}, true},
		{"dropped", Result{Confirmed: []*recordv1.Record{a}}, false},
		{"duplicate", Result{Confirmed: batch, Remaining: []*recordv1.Record{a}}, false},
		{"extra record", Result{Confirmed: append([]*recordv1.Record{a, b}, newScalar(3, "c"))}, false},
		{"nothing", Result{}, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := conserved(batch, tc.result); got != tc.want {
				t.Fatalf("conserved = %v, want %v", got, tc.want)
			}
		})
	}
}
