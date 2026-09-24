package server

import (
	"testing"
	"time"

	"google.golang.org/grpc"
)

const testGrace = 2 * time.Second

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
