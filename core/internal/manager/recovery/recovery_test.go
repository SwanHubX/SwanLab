package recovery

import (
	"context"
	"sync"
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

func TestPumpStopBeforeStartAndConcurrentLifecycle(t *testing.T) {
	// Stop -> Start -> Stop：先 Stop 后，再 Start 不得执行 DrainRecovery。
	for range 16 {
		src := &fakeSource{}
		var degraded atomic.Bool
		degraded.Store(true)
		p := New(Options{
			Source:   src,
			Degraded: &degraded,
			Interval: 2 * time.Millisecond,
		})
		p.Stop()
		p.Start()
		p.Stop()
		p.Stop()
		if got := p.lifecycle.Load(); got != lifecycleStopped {
			t.Fatalf("lifecycle = %d after Stop -> Start, want lifecycleStopped (%d)", got, lifecycleStopped)
		}
		if got := src.calls.Load(); got != 0 {
			t.Fatalf("DrainRecovery calls = %d after Stop -> Start, want 0", got)
		}
	}

	// 并发 Start / Stop 在超时预算内稳定退出，且未启用 Pump 安全空调用。
	for range 24 {
		src := &fakeSource{}
		var degraded atomic.Bool
		p := New(Options{
			Source:   src,
			Degraded: &degraded,
			Interval: 2 * time.Millisecond,
		})
		var wg sync.WaitGroup
		wg.Add(4)
		go func() { defer wg.Done(); p.Start() }()
		go func() { defer wg.Done(); p.Start() }()
		go func() { defer wg.Done(); p.Stop() }()
		go func() { defer wg.Done(); p.Stop() }()

		done := make(chan struct{})
		go func() {
			wg.Wait()
			close(done)
		}()
		select {
		case <-done:
		case <-time.After(2 * time.Second):
			t.Fatal("concurrent Pump Start/Stop timed out")
		}
		p.Stop()
	}

	disabled := New(Options{})
	disabled.Start()
	disabled.Stop()
}
