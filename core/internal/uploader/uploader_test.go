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
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
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

// successUpload 是恒成功的 Upload 实现。
func successUpload(_ context.Context, records []*recordv1.Record) Result {
	return Result{Confirmed: records}
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
	if _, err := New(Options{
		In:      make(chan *recordv1.Record),
		Budget:  ledger,
		Senders: map[record.Type]Sender{record.TypeScalar: nil},
	}); err == nil {
		t.Fatal("nil Senders should fail")
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
	if batch := recvBatch(t, calls); len(batch) != 2 {
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
	if batch := recvBatch(t, calls); len(batch) != 2 {
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
	flushDone := make(chan error, 1)
	go func() { flushDone <- f.u.FlushSchedulingForTest() }()
	waitSignal(t, started, "slow sender start")
	if f.u.Drained() {
		t.Fatal("drained must stay false while a batch is in flight")
	}

	close(f.in)
	if f.u.Drained() {
		t.Fatal("drained must stay false while the sender is in flight")
	}
	close(release)
	if err := <-flushDone; err != nil {
		t.Fatalf("flush: %v", err)
	}
	waitDone(t, f.u)
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

func TestColumnCompletesBeforeScalar(t *testing.T) {
	started := make(chan struct{}, 1)
	release := make(chan struct{})
	order := make(chan record.Type, 2)
	column := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		started <- struct{}{}
		<-release
		order <- record.TypeColumn
		return Result{Confirmed: records}
	})
	scalar := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		order <- record.TypeScalar
		return Result{Confirmed: records}
	})
	f := newUplFixture(t, map[record.Type]Sender{
		record.TypeColumn: column,
		record.TypeScalar: scalar,
	}, 100, time.Hour)

	f.feed(t, newScalar(2, "loss"), newColumn(1, "loss"))
	flushDone := make(chan error, 1)
	go func() { flushDone <- f.u.FlushSchedulingForTest() }()
	waitSignal(t, started, "column upload")
	select {
	case typ := <-order:
		t.Fatalf("%s ran while column was in flight", typ)
	default:
	}
	close(release)
	if got := <-order; got != record.TypeColumn {
		t.Fatalf("first completed sender = %s, want column", got)
	}
	if got := <-order; got != record.TypeScalar {
		t.Fatalf("second completed sender = %s, want scalar", got)
	}
	if err := <-flushDone; err != nil {
		t.Fatalf("flush: %v", err)
	}
	close(f.in)
	waitDone(t, f.u)
}

func TestFailureRetainsRemainingAndLaterGroups(t *testing.T) {
	scalar := &scriptedSender{scripts: []func([]*recordv1.Record) Result{
		func(records []*recordv1.Record) Result {
			return Result{Confirmed: records[:1], Remaining: records[1:]}
		},
		func(records []*recordv1.Record) Result { return Result{Confirmed: records} },
	}}
	logCalls := 0
	logSender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		logCalls++
		return Result{Confirmed: records}
	})
	f := newUplFixture(t, map[record.Type]Sender{
		record.TypeColumn: senderFunc(successUpload),
		record.TypeScalar: scalar,
		record.TypeLog:    logSender,
	}, 100, time.Hour)

	f.feed(t, newColumn(1, "loss"), newScalar(2, "a"), newScalar(3, "b"), newLog(4, "line"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	if logCalls != 0 {
		t.Fatal("later log group must wait for the failed scalar group")
	}
	stats := f.u.Stats()
	if stats.Pending != 2 || stats.TotalRecords != 4 || stats.UploadedRecords != 2 {
		t.Fatalf("stats = %+v, want pending=2 total=4 uploaded=2", stats)
	}
	if used, _ := f.ledger.Usage(); used != 2 {
		t.Fatalf("budget usage = %d, want 2", used)
	}

	f.u.RequestFinish()
	close(f.in)
	waitDone(t, f.u)
	if logCalls != 1 {
		t.Fatalf("log calls = %d, want 1 after scalar retry", logCalls)
	}
	if got := f.tracker.Snapshot().UploadedRecords; got != 4 {
		t.Fatalf("uploaded records = %d, want 4", got)
	}
	if used, _ := f.ledger.Usage(); used != 0 {
		t.Fatalf("budget usage = %d, want 0", used)
	}
}

func TestInvalidResultRetriesWholeBatch(t *testing.T) {
	sender := &scriptedSender{scripts: []func([]*recordv1.Record) Result{
		func(records []*recordv1.Record) Result { return Result{Confirmed: records[:1]} },
		func(records []*recordv1.Record) Result { return Result{Confirmed: records} },
	}}
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 100, time.Hour)
	f.feed(t, newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	if got := f.u.Stats().Pending; got != 3 {
		t.Fatalf("pending = %d, want whole batch", got)
	}
	if failures := f.u.LastFailures(); len(failures) != 1 || failures[0].Message != "sender result does not account for every record" {
		t.Fatalf("failures = %+v, want conservation failure", failures)
	}

	f.u.RequestFinish()
	close(f.in)
	waitDone(t, f.u)
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

	attempt := 0
	confirmed := make(map[*recordv1.Record]int, total)
	sender := senderFunc(func(_ context.Context, batch []*recordv1.Record) Result {
		attempt++
		first := attempt == 1
		var done, retry []*recordv1.Record
		for i, rec := range batch {
			if first && i%2 == 1 {
				retry = append(retry, rec)
				continue
			}
			done = append(done, rec)
		}
		for _, rec := range done {
			confirmed[rec]++
		}
		return Result{Confirmed: done, Remaining: retry}
	})
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 7, time.Hour)
	f.feed(t, recs...)
	close(f.in)
	waitDone(t, f.u)

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
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 100, time.Hour)
	f.u.retry = retryNoJitter(time.Hour)
	f.feed(t, newScalar(1, "a"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	<-attempts

	f.u.RequestFinish()
	select {
	case <-attempts:
	case <-time.After(2 * time.Second):
		t.Fatal("RequestFinish did not interrupt the backoff window")
	}
	close(f.in)
	waitDone(t, f.u)
}

func TestRequestFinishClearsBackoffOnlyOnce(t *testing.T) {
	attempts := make(chan int, 4)
	var count atomic.Int32
	sender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		attempts <- int(count.Add(1))
		return Result{Remaining: records}
	})
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 100, time.Hour)
	f.u.retry = retryNoJitter(time.Hour)
	f.feed(t, newScalar(1, "a"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	<-attempts

	f.u.RequestFinish()
	select {
	case <-attempts:
	case <-time.After(2 * time.Second):
		t.Fatal("first RequestFinish did not interrupt backoff")
	}
	f.u.RequestFinish()
	select {
	case attempt := <-attempts:
		t.Fatalf("attempt %d ran without a new backoff after repeated RequestFinish", attempt)
	case <-time.After(40 * time.Millisecond):
	}
	if got := count.Load(); got != 2 {
		t.Fatalf("attempts = %d, want bounded count 2", got)
	}
}

func TestSuspendResume(t *testing.T) {
	calls, sender := recordingSender(1)
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 100, time.Hour)
	f.u.Suspend()
	f.feed(t, newScalar(1, "a"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush while suspended: %v", err)
	}
	select {
	case <-calls:
		t.Fatal("sender ran while suspended")
	case <-time.After(40 * time.Millisecond):
	}

	f.u.Resume()
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush after resume: %v", err)
	}
	if batch := recvBatch(t, calls); len(batch) != 1 {
		t.Fatalf("batch len = %d, want 1", len(batch))
	}
	close(f.in)
	waitDone(t, f.u)
}

func TestStatsAndSnapshot(t *testing.T) {
	calls, sender := recordingSender(1)
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 3, time.Hour)
	f.feed(t, newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c"))
	recvBatch(t, calls)
	close(f.in)
	waitDone(t, f.u)

	stats := f.u.Stats()
	if stats.Buffered != 0 || stats.Pending != 0 || stats.TotalRecords != 3 || stats.UploadedRecords != 3 {
		t.Fatalf("stats = %+v, want empty queues and 3/3 records", stats)
	}
	f.tracker.SetState(operationv1.CoreState_CORE_STATE_RUNNING)
	f.tracker.AddFile("f1", "/tmp/a", 100)
	f.tracker.UpdateFileProgress("f1", "part", 40)
	snapshot := f.u.Snapshot()
	if snapshot.State != operationv1.CoreState_CORE_STATE_RUNNING ||
		snapshot.TotalNumber != 1 || snapshot.TotalSize != 100 || snapshot.UploadedSize != 40 ||
		snapshot.TotalRecords != 3 || snapshot.UploadedRecords != 3 || len(snapshot.Files) != 1 {
		t.Fatalf("snapshot = %+v, want complete OperationStats", snapshot)
	}
}

func TestTerminalFailureWarns(t *testing.T) {
	sender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		return Result{
			Confirmed: records,
			Failures:  []Failure{{Message: "upload denied", Records: records}},
		}
	})
	f := newUplFixture(t, map[record.Type]Sender{record.TypeScalar: sender}, 100, time.Hour)
	var warnings []string
	f.u.warnFailure = func(message string) { warnings = append(warnings, message) }
	f.feed(t, newScalar(1, "a"))
	if err := f.u.FlushSchedulingForTest(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	if len(warnings) != 1 || warnings[0] != "upload denied" {
		t.Fatalf("warnings = %v, want terminal failure warning", warnings)
	}
	close(f.in)
	waitDone(t, f.u)
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
