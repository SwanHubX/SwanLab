package server

import (
	"context"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/stats"
)

type fakeTimer struct {
	mu      sync.Mutex
	stopped bool
	fn      func()
}

func (f *fakeTimer) Stop() bool {
	f.mu.Lock()
	defer f.mu.Unlock()
	wasActive := !f.stopped
	f.stopped = true
	return wasActive
}

func (f *fakeTimer) isStopped() bool {
	f.mu.Lock()
	defer f.mu.Unlock()
	return f.stopped
}

func (f *fakeTimer) fire() {
	f.fn()
}

type fakeClockHarness struct {
	mu     sync.Mutex
	timers []*fakeTimer
}

func (h *fakeClockHarness) after(_ time.Duration, fn func()) Timer {
	t := &fakeTimer{fn: fn}
	h.mu.Lock()
	h.timers = append(h.timers, t)
	h.mu.Unlock()
	return t
}

func (h *fakeClockHarness) count() int {
	h.mu.Lock()
	defer h.mu.Unlock()
	return len(h.timers)
}

func (h *fakeClockHarness) last() *fakeTimer {
	h.mu.Lock()
	defer h.mu.Unlock()
	return h.timers[len(h.timers)-1]
}

func TestIdleMonitorStartsAtZeroAndFiresOnce(t *testing.T) {
	harness := &fakeClockHarness{}
	var fired atomic.Int32
	m := NewIdleMonitor(2*time.Second, func() { fired.Add(1) }, WithAfterFunc(harness.after))

	m.Start()
	if harness.count() != 1 {
		t.Fatalf("timers created = %d, want 1", harness.count())
	}

	harness.last().fire()
	// 重复触发同一批回调不会二次调用 onIdle。
	harness.last().fire()
	if fired.Load() != 1 {
		t.Fatalf("onIdle calls = %d, want 1", fired.Load())
	}
}

func TestIdleMonitorCancelledByConnAndRearmedOnDisconnect(t *testing.T) {
	harness := &fakeClockHarness{}
	var fired atomic.Int32
	m := NewIdleMonitor(2*time.Second, func() { fired.Add(1) }, WithAfterFunc(harness.after))

	m.Start()
	first := harness.last()

	// 新连接进入时取消首个空闲定时器。
	m.OnConnBegin()
	if !first.isStopped() {
		t.Fatal("first timer must be stopped on new connection")
	}
	if m.Active() != 1 {
		t.Fatalf("active = %d, want 1", m.Active())
	}

	// 即使旧定时器并发到期，也不应误触发关闭。
	first.fire()
	if fired.Load() != 0 {
		t.Fatalf("stale timer fired onIdle (%d times)", fired.Load())
	}

	// 还有第二条连接时，断开第一条不应启动空闲计时。
	m.OnConnBegin()
	m.OnConnEnd()
	if harness.count() != 1 {
		t.Fatalf("timers = %d, must not arm while active=1", harness.count())
	}

	// 最后一条连接断开后重新计时，超时触发关闭。
	m.OnConnEnd()
	if harness.count() != 2 {
		t.Fatalf("timers = %d, want 2 after last conn closes", harness.count())
	}
	harness.last().fire()
	if fired.Load() != 1 {
		t.Fatalf("onIdle calls = %d, want 1", fired.Load())
	}
}

func TestIdleMonitorDisabledWhenZeroTimeout(t *testing.T) {
	harness := &fakeClockHarness{}
	m := NewIdleMonitor(0, func() {}, WithAfterFunc(harness.after))
	m.Start()
	m.OnConnBegin()
	m.OnConnEnd()
	if harness.count() != 0 {
		t.Fatalf("zero timeout must not schedule timer, got %d", harness.count())
	}
}

func TestConnTrackerDeduplicatesConnEnd(t *testing.T) {
	harness := &fakeClockHarness{}
	var fired atomic.Int32
	m := NewIdleMonitor(time.Second, func() { fired.Add(1) }, WithAfterFunc(harness.after))
	tracker := NewConnTracker()
	tracker.SetMonitor(m)
	m.Start()

	connCtx := tracker.TagConn(context.Background(), &stats.ConnTagInfo{})
	tracker.HandleConn(connCtx, &stats.ConnBegin{})
	if tracker.Active() != 1 || m.Active() != 1 {
		t.Fatalf("active tracker=%d monitor=%d, want 1", tracker.Active(), m.Active())
	}

	tracker.HandleConn(connCtx, &stats.ConnEnd{})
	tracker.HandleConn(connCtx, &stats.ConnEnd{})
	if tracker.Active() != 0 || m.Active() != 0 {
		t.Fatalf("duplicate ConnEnd skewed active count: tracker=%d monitor=%d", tracker.Active(), m.Active())
	}

	harness.last().fire()
	if fired.Load() != 1 {
		t.Fatalf("onIdle calls = %d, want 1", fired.Load())
	}
}

func TestControllerStopsIdleMonitorOnExternalShutdown(t *testing.T) {
	harness := &fakeClockHarness{}
	g := grpc.NewServer()
	ctrl := NewController(g, testGrace)
	m := NewIdleMonitor(time.Second, func() { ctrl.Shutdown("idle-timeout") }, WithAfterFunc(harness.after))
	ctrl.SetIdleMonitor(m)
	m.Start()

	timer := harness.last()
	ctrl.Shutdown("signal")
	<-ctrl.Done()

	if !timer.isStopped() {
		t.Fatal("Controller.Shutdown must stop idle timer")
	}
}
