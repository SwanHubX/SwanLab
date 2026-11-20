package recovery

import (
	"context"
	"sync/atomic"
	"testing"
	"time"

	"go.uber.org/goleak"

	"github.com/swanhubx/swanlab/core/internal/buffer"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

type fakeSource struct {
	calls     atomic.Int32
	remaining atomic.Int32
	err       error
}

func (f *fakeSource) DrainRecovery(context.Context) (int, error) {
	n := f.calls.Add(1)
	if f.err != nil {
		return int(f.remaining.Load()), f.err
	}
	if n >= 3 {
		f.remaining.Store(0)
	}
	return int(f.remaining.Load()), nil
}

func TestPumpSkipsUntilDegradedAndStopsCleanly(t *testing.T) {
	src := &fakeSource{}
	var degraded atomic.Bool
	p := New(Options{
		Source:   src,
		Degraded: &degraded,
		Interval: 5 * time.Millisecond,
	})
	p.Start()
	time.Sleep(25 * time.Millisecond)
	if got := src.calls.Load(); got != 0 {
		t.Fatalf("calls before degraded = %d, want 0", got)
	}

	degraded.Store(true)
	deadline := time.Now().Add(time.Second)
	for src.calls.Load() == 0 && time.Now().Before(deadline) {
		time.Sleep(5 * time.Millisecond)
	}
	if src.calls.Load() == 0 {
		t.Fatal("expected recovery drain call once degraded")
	}

	p.Stop()
	p.Stop() // 幂等重复调用
}

func TestPumpExitsOnClosedError(t *testing.T) {
	src := &fakeSource{err: buffer.ErrClosed}
	var degraded atomic.Bool
	degraded.Store(true)
	p := New(Options{
		Source:   src,
		Degraded: &degraded,
		Interval: 5 * time.Millisecond,
	})
	p.Start()
	select {
	case <-p.done:
	case <-time.After(time.Second):
		t.Fatal("pump did not exit after ErrClosed")
	}
	p.Stop()
}

func TestPumpDrainPollsUntilEmptyOrTimeout(t *testing.T) {
	src := &fakeSource{}
	src.remaining.Store(2)
	p := New(Options{
		Source:       src,
		PollInterval: 5 * time.Millisecond,
	})
	p.Drain(context.Background())
	if got := src.calls.Load(); got != 3 {
		t.Fatalf("Drain calls = %d, want 3", got)
	}

	// nil receiver safety
	var nilPump *Pump
	nilPump.Start()
	nilPump.Drain(context.Background())
	nilPump.Stop()

	// ctx timeout safety when remaining never drops to 0
	stuck := &fakeSource{}
	stuck.remaining.Store(5)
	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
	defer cancel()
	New(Options{Source: stuck, PollInterval: 5 * time.Millisecond}).Drain(ctx)
}
