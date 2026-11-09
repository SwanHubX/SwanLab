package uploader

import (
	"context"
	"sync/atomic"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

func TestDispatchCallbacksReportRounds(t *testing.T) {
	ledger := budget.New(budget.Limits{})
	tracker := NewTracker()
	in := make(chan *recordv1.Record, 8)
	rounds := make(chan int, 8)
	success := make(chan struct{}, 8)

	var fail atomic.Bool
	fail.Store(true)
	sender := senderFunc(func(_ context.Context, records []*recordv1.Record) Result {
		if fail.Load() {
			return Result{Remaining: records}
		}
		return Result{Confirmed: records}
	})

	u, err := New(Options{
		In:       in,
		Budget:   ledger,
		Senders:  map[record.Type]Sender{record.TypeScalar: sender},
		Batch:    1,
		Interval: time.Hour,
		Retry:    retryNoJitter(time.Millisecond),
		Tracker:  tracker,
		Warn:     quietThrottle(),
		OnRetryableFailure: func(roundsDone int) {
			rounds <- roundsDone
		},
		OnDispatchSuccess: func() { success <- struct{}{} },
	})
	if err != nil {
		t.Fatalf("New: %v", err)
	}
	t.Cleanup(func() { _ = u.Close() })

	rec := newScalar(1, "a")
	tracker.AddTotalRecords(1)
	if err := ledger.Acquire(context.Background(), 1, record.Size(rec)); err != nil {
		t.Fatalf("acquire: %v", err)
	}
	in <- rec

	for want := 1; want <= 3; want++ {
		select {
		case got := <-rounds:
			if got != want {
				t.Fatalf("round = %d, want %d", got, want)
			}
		case <-time.After(2 * time.Second):
			t.Fatalf("timeout waiting for failure round %d", want)
		}
	}

	fail.Store(false)
	select {
	case <-success:
	case <-time.After(2 * time.Second):
		t.Fatal("timeout waiting for dispatch success")
	}
}
