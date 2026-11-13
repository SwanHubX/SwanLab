package cron

import (
	"context"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"go.uber.org/goleak"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

type fakeClock struct {
	mu      sync.Mutex
	chans   map[time.Duration]chan time.Time
	changed chan struct{}
}

func newFakeClock() *fakeClock {
	return &fakeClock{
		chans:   make(map[time.Duration]chan time.Time),
		changed: make(chan struct{}, 1),
	}
}

func (c *fakeClock) Ticker(d time.Duration) (<-chan time.Time, func()) {
	ch := make(chan time.Time, 8)
	c.mu.Lock()
	c.chans[d] = ch
	c.mu.Unlock()
	select {
	case c.changed <- struct{}{}:
	default:
	}
	return ch, func() {}
}

// ticks 等待指定周期的 ticker 注册，避免依赖协程启动顺序。
func (c *fakeClock) ticks(t *testing.T, d time.Duration) chan time.Time {
	t.Helper()
	deadline := time.After(2 * time.Second)
	for {
		c.mu.Lock()
		ch, ok := c.chans[d]
		c.mu.Unlock()
		if ok {
			return ch
		}
		select {
		case <-c.changed:
		case <-deadline:
			t.Fatalf("ticker for %v was not created", d)
		}
	}
}

func waitRun(t *testing.T, runs <-chan string) string {
	t.Helper()
	select {
	case job := <-runs:
		return job
	case <-time.After(2 * time.Second):
		t.Fatal("timeout waiting for job run")
		return ""
	}
}

// 同一任务阻塞期间的多个 tick 合并为一次执行，不追帧。
func TestSchedulerCoalescesBacklog(t *testing.T) {
	clock := newFakeClock()
	started := make(chan struct{})
	release := make(chan struct{})
	runs := make(chan string, 4)
	var calls atomic.Int32

	s := New(clock)
	s.Add(Job{Interval: time.Second, Run: func(context.Context) {
		if calls.Add(1) == 1 {
			close(started)
			<-release
		}
		runs <- "job"
	}})
	s.Start(context.Background())
	ticks := clock.ticks(t, time.Second)

	ticks <- time.Time{}
	select {
	case <-started:
	case <-time.After(2 * time.Second):
		t.Fatal("job did not start")
	}
	// 任务阻塞期间到达的 3 个 tick 必须被合并。
	for range 3 {
		ticks <- time.Time{}
	}
	close(release)

	waitRun(t, runs)
	waitRun(t, runs)
	if got := calls.Load(); got != 2 {
		t.Fatalf("runs = %d, want 2 after coalescing", got)
	}
	s.Stop()
	s.Stop()
}

func TestSchedulerRunsJobsIndependently(t *testing.T) {
	clock := newFakeClock()
	runs := make(chan string, 4)

	s := New(clock)
	s.Add(Job{Interval: time.Second, Run: func(context.Context) { runs <- "first" }})
	s.Add(Job{Interval: time.Minute, Run: func(context.Context) { runs <- "second" }})
	s.Start(context.Background())
	fast := clock.ticks(t, time.Second)
	slow := clock.ticks(t, time.Minute)

	fast <- time.Time{}
	if job := waitRun(t, runs); job != "first" {
		t.Fatalf("job = %q, want first", job)
	}
	select {
	case job := <-runs:
		t.Fatalf("second job ran without a tick: %q", job)
	default:
	}

	slow <- time.Time{}
	if job := waitRun(t, runs); job != "second" {
		t.Fatalf("job = %q, want second", job)
	}
	s.Stop()
}

func TestSchedulerValidatesJobs(t *testing.T) {
	s := New(newFakeClock())
	for name, job := range map[string]Job{
		"missing interval": {Run: func(context.Context) {}},
		"missing run":      {Interval: time.Second},
	} {
		t.Run(name, func(t *testing.T) {
			defer func() {
				if recover() == nil {
					t.Fatal("expected panic")
				}
			}()
			s.Add(job)
		})
	}

	s.Add(Job{Interval: time.Second, Run: func(context.Context) {}})
	s.Start(context.Background())
	defer s.Stop()
	defer func() {
		if recover() == nil {
			t.Fatal("Add after Start should panic")
		}
	}()
	s.Add(Job{Interval: time.Second, Run: func(context.Context) {}})
}
