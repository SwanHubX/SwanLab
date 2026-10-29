package uploader

import (
	"context"
	"errors"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// recordingSender 返回一个把每次调用批次写入通道的 sender。
func recordingSender(size int) (chan []*recordv1.Record, Sender) {
	calls := make(chan []*recordv1.Record, size)
	sender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		calls <- records
		return Result{Confirmed: records}
	})
	return calls, sender
}

func recvBatch(t *testing.T, calls <-chan []*recordv1.Record) []*recordv1.Record {
	t.Helper()
	select {
	case batch := <-calls:
		return batch
	case <-time.After(2 * time.Second):
		t.Fatal("timeout waiting for sender call")
		return nil
	}
}

type uplFixture struct {
	u       *Uploader
	in      chan *recordv1.Record
	ledger  *budget.Budget
	tracker *Tracker
}

func newUplFixture(t *testing.T, senders map[record.Type]Sender, batch int, interval time.Duration) *uplFixture {
	t.Helper()
	ledger := budget.New(budget.Limits{})
	tracker := NewTracker()
	in := make(chan *recordv1.Record, 256)
	u, err := New(Options{
		In:       in,
		Budget:   ledger,
		Senders:  senders,
		Batch:    batch,
		Interval: interval,
		Retry:    retryNoJitter(time.Millisecond),
		Tracker:  tracker,
		Warn:     quietThrottle(),
	})
	if err != nil {
		t.Fatalf("new uploader: %v", err)
	}
	t.Cleanup(func() { _ = u.Close() })
	return &uplFixture{u: u, in: in, ledger: ledger, tracker: tracker}
}

func (f *uplFixture) feed(t *testing.T, records ...*recordv1.Record) {
	t.Helper()
	for _, rec := range records {
		if err := f.ledger.Acquire(context.Background(), 1, record.Size(rec)); err != nil {
			t.Fatalf("acquire: %v", err)
		}
		f.in <- rec
	}
}

func waitDone(t *testing.T, u *Uploader) {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), 2*time.Second)
	defer cancel()
	if err := u.Wait(ctx); err != nil {
		t.Fatalf("wait for drain: %v", err)
	}
}

func TestNewRequiresCoreOptions(t *testing.T) {
	ledger := budget.New(budget.Limits{})
	if _, err := New(Options{}); err == nil {
		t.Fatal("missing In should fail")
	}
	if _, err := New(Options{In: make(chan *recordv1.Record)}); err == nil {
		t.Fatal("missing Budget should fail")
	}
	if _, err := New(Options{In: make(chan *recordv1.Record), Budget: ledger}); err == nil {
		t.Fatal("missing Senders should fail")
	}
}

func TestBatchSizeTriggersFlush(t *testing.T) {
	calls, sender := recordingSender(2)
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 3, time.Hour)

	recs := []*recordv1.Record{newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c")}
	f.feed(t, recs...)

	batch := recvBatch(t, calls)
	if len(batch) != 3 || batch[0] != recs[0] || batch[1] != recs[1] || batch[2] != recs[2] {
		t.Fatalf("batch = %v, want all three in order", batch)
	}

	close(f.in)
	waitDone(t, f.u)
	if !f.u.Drained() {
		t.Fatal("Drained() = false after drain")
	}
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0", used)
	}
}

func TestIntervalTriggersFlush(t *testing.T) {
	calls, sender := recordingSender(2)
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 100, 50*time.Millisecond)

	f.feed(t, newScalar(1, "a"), newScalar(2, "b"))

	batch := recvBatch(t, calls)
	if len(batch) != 2 {
		t.Fatalf("batch len = %d, want 2", len(batch))
	}

	close(f.in)
	waitDone(t, f.u)
}

func TestFlushSchedulingForTest(t *testing.T) {
	calls, sender := recordingSender(2)
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 100, time.Hour)

	f.feed(t, newScalar(1, "a"), newScalar(2, "b"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}

	batch := recvBatch(t, calls)
	if len(batch) != 2 {
		t.Fatalf("batch len = %d, want 2", len(batch))
	}

	close(f.in)
	waitDone(t, f.u)
	if err := f.u.FlushSchedulingForTest(); !errors.Is(err, ErrClosed) {
		t.Fatalf("flush after close = %v, want ErrClosed", err)
	}
}

func TestDrainedNotReportedEarly(t *testing.T) {
	started := make(chan struct{}, 1)
	release := make(chan struct{})
	slow := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		started <- struct{}{}
		<-release
		return Result{Confirmed: records}
	})
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: slow}, 100, time.Hour)

	f.feed(t, newScalar(1, "a"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	waitSignal(t, started, "slow sender start")

	if f.u.Drained() {
		t.Fatal("drained must stay false while a batch is in flight")
	}

	close(f.in)
	if f.u.Drained() {
		t.Fatal("drained must stay false while the entry is still open/in flight")
	}

	close(release)
	waitDone(t, f.u)
	if !f.u.Drained() {
		t.Fatal("drained should be true after the in-flight batch completes")
	}
}

func TestOnDrainedCalledOnce(t *testing.T) {
	ledger := budget.New(budget.Limits{})
	in := make(chan *recordv1.Record)
	var calls atomic.Int32
	u, err := New(Options{
		In:       in,
		Budget:   ledger,
		Senders:  map[record.Type]Sender{record.TypeScalar: senderFunc(successUpload)},
		Batch:    10,
		Interval: time.Hour,
		Warn:     quietThrottle(),
		OnDrained: func() {
			calls.Add(1)
		},
	})
	if err != nil {
		t.Fatalf("new uploader: %v", err)
	}
	t.Cleanup(func() { _ = u.Close() })

	close(in)
	waitDone(t, u)
	if got := calls.Load(); got != 1 {
		t.Fatalf("OnDrained calls = %d, want 1", got)
	}
}

// successUpload 是恒成功的 Upload 实现。
func successUpload(_ context.Context, records []*recordv1.Record) Result {
	return Result{Confirmed: records}
}

func TestRequestFinishInterruptsBackoff(t *testing.T) {
	attempts := make(chan int, 4)
	var count atomic.Int32
	sender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		attempt := int(count.Add(1))
		attempts <- attempt
		if attempt == 1 {
			return Result{Remaining: records}
		}
		return Result{Confirmed: records}
	})

	ledger := budget.New(budget.Limits{})
	in := make(chan *recordv1.Record, 16)
	u, err := New(Options{
		In:       in,
		Budget:   ledger,
		Senders:  map[record.Type]Sender{record.TypeScalar: sender},
		Batch:    100,
		Interval: time.Hour,
		Retry:    retryNoJitter(time.Hour),
		Warn:     quietThrottle(),
	})
	if err != nil {
		t.Fatalf("new uploader: %v", err)
	}
	t.Cleanup(func() { _ = u.Close() })

	rec := newScalar(1, "a")
	if err := ledger.Acquire(context.Background(), 1, record.Size(rec)); err != nil {
		t.Fatalf("acquire: %v", err)
	}
	in <- rec
	if err := u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	select {
	case <-attempts:
	case <-time.After(2 * time.Second):
		t.Fatal("timeout waiting for first attempt")
	}

	u.RequestFinish()
	select {
	case <-attempts:
	case <-time.After(2 * time.Second):
		t.Fatal("RequestFinish did not interrupt the backoff window")
	}

	close(in)
	waitDone(t, u)
}

func TestCloseIdempotent(t *testing.T) {
	t.Run("without drain", func(t *testing.T) {
		f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: senderFunc(successUpload)}, 10, time.Hour)
		if err := f.u.Close(); !errors.Is(err, ErrNotDrained) {
			t.Fatalf("first close = %v, want ErrNotDrained", err)
		}
		if err := f.u.Close(); !errors.Is(err, ErrNotDrained) {
			t.Fatalf("second close = %v, want ErrNotDrained", err)
		}
	})

	t.Run("after drain", func(t *testing.T) {
		f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: senderFunc(successUpload)}, 10, time.Hour)
		close(f.in)
		waitDone(t, f.u)
		if err := f.u.Close(); err != nil {
			t.Fatalf("close after drain = %v, want nil", err)
		}
		if err := f.u.Close(); err != nil {
			t.Fatalf("second close = %v, want nil", err)
		}
	})
}

func TestWaitHonorsContext(t *testing.T) {
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: senderFunc(successUpload)}, 10, time.Hour)

	ctx, cancel := context.WithTimeout(context.Background(), 30*time.Millisecond)
	defer cancel()
	if err := f.u.Wait(ctx); !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("wait = %v, want DeadlineExceeded", err)
	}
}

func TestStats(t *testing.T) {
	calls, sender := recordingSender(1)
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 3, time.Hour)

	f.feed(t, newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c"))
	recvBatch(t, calls)
	close(f.in)
	waitDone(t, f.u)

	stats := f.u.Stats()
	if stats.Buffered != 0 || stats.TotalRecords != 3 || stats.UploadedRecords != 3 {
		t.Fatalf("stats = %+v, want buffered=0 total=3 uploaded=3", stats)
	}
	if len(stats.Tracks) != 1 {
		t.Fatalf("tracks = %d, want 1", len(stats.Tracks))
	}
	track := stats.Tracks[0]
	if track.Type != record.TypeScalar || track.Queue != 0 || track.Confirmed != 3 {
		t.Fatalf("track stats = %+v", track)
	}
}

func TestUnregisteredTypeTerminates(t *testing.T) {
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: senderFunc(successUpload)}, 100, time.Hour)

	var mu sync.Mutex
	var warned []record.Type
	f.u.warnUnsupported = func(typ record.Type, _ int) {
		mu.Lock()
		warned = append(warned, typ)
		mu.Unlock()
	}

	f.feed(t, newLog(1, "line-1"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	f.feed(t, newLog(2, "line-2"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}

	close(f.in)
	waitDone(t, f.u)

	mu.Lock()
	defer mu.Unlock()
	if len(warned) != 1 || warned[0] != record.TypeLog {
		t.Fatalf("warnings = %v, want exactly one log warning", warned)
	}
	stats := f.u.Stats()
	if stats.TotalRecords != 2 || stats.UploadedRecords != 2 {
		t.Fatalf("stats = %+v, want both log records terminal", stats)
	}
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0", used)
	}
}
