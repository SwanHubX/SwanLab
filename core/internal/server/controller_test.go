package server

import (
	"context"
	"errors"
	"sync/atomic"
	"testing"
	"time"

	"go.uber.org/goleak"
	"google.golang.org/grpc"
)

const testGrace = 2 * time.Second

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

func TestControllerRunsShutdownHook(t *testing.T) {
	g := grpc.NewServer()
	ctrl := NewController(g, testGrace)
	var calls atomic.Int32
	ctrl.SetOnShutdown(func(ctx context.Context) error {
		calls.Add(1)
		if _, ok := ctx.Deadline(); !ok {
			return errors.New("shutdown hook context must carry the shared deadline")
		}
		return nil
	})
	ctrl.Shutdown("hook")
	select {
	case <-ctrl.Done():
	case <-time.After(3 * time.Second):
		t.Fatal("controller Done not closed after hook")
	}
	if got := calls.Load(); got != 1 {
		t.Fatalf("shutdown hook called %d times, want 1", got)
	}
}

// hook 与 GracefulStop 必须共享同一个 grace deadline：hook 等满预算后
// 控制器仍要进入 Stop 兜底，总耗时只受一个 grace 预算约束。
func TestControllerShutdownSharesGraceWithHook(t *testing.T) {
	const grace = 300 * time.Millisecond
	g := grpc.NewServer()
	ctrl := NewController(g, grace)
	entered := make(chan struct{})
	ctrl.SetOnShutdown(func(ctx context.Context) error {
		close(entered)
		<-ctx.Done()
		return ctx.Err()
	})

	start := time.Now()
	ctrl.Shutdown("timeout")
	<-entered
	select {
	case <-ctrl.Done():
	case <-time.After(3 * time.Second):
		t.Fatal("controller Done not closed after timed-out hook")
	}
	elapsed := time.Since(start)
	if elapsed < grace/2 {
		t.Fatalf("hook returned too early, elapsed=%s", elapsed)
	}
	if elapsed > 3*grace {
		t.Fatalf("shutdown hook and GracefulStop must share one deadline, elapsed=%s grace=%s", elapsed, grace)
	}
}

func TestControllerShutdownIdempotent(t *testing.T) {
	g := grpc.NewServer()
	ctrl := NewController(g, testGrace)
	ctrl.Shutdown("first")
	ctrl.Shutdown("second")
	select {
	case <-ctrl.Done():
	case <-time.After(2 * time.Second):
		t.Fatal("controller Done not closed after Shutdown")
	}
	<-ctrl.Done()
	if st := ctrl.Lifecycle().Get(); st != StateClosed {
		t.Fatalf("lifecycle state after shutdown = %v, want StateClosed", st)
	}
}

func TestLifecycleStateTransitions(t *testing.T) {
	lc := newLifecycle()
	if lc.Get() != StateNotReady {
		t.Fatalf("initial state = %v, want StateNotReady", lc.Get())
	}
	if !lc.Spinup() || lc.Get() != StateReady {
		t.Fatalf("spinup state = %v, want StateReady", lc.Get())
	}
	if !lc.Spinup() {
		t.Fatal("repeated Spinup in READY must succeed")
	}
	lc.BeginStopping()
	if lc.Get() != StateStopping || lc.Spinup() {
		t.Fatalf("stopping state = %v, Spinup must be rejected", lc.Get())
	}
	lc.Close()
	lc.Close()
	if lc.Get() != StateClosed {
		t.Fatalf("closed state = %v, want StateClosed", lc.Get())
	}
}
