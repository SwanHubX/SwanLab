package manager

import (
	"context"
	"errors"
	"fmt"
	"sync"
	"sync/atomic"
	"testing"

	"github.com/google/uuid"
)

func mustState(t *testing.T, m *Manager, handle string, want RunState) {
	t.Helper()
	got, err := m.State(handle)
	if err != nil {
		t.Fatalf("State(%q): %v", handle, err)
	}
	if got != want {
		t.Fatalf("State(%q) = %s, want %s", handle, got, want)
	}
}

func TestManagerLifecycle(t *testing.T) {
	m := New()
	ctx := context.Background()
	handle, err := m.Start()
	if err != nil || handle == "" {
		t.Fatalf("Start: handle=%q err=%v", handle, err)
	}
	mustState(t, m, handle, RunStateRunning)

	if err := m.Confirm(ctx, handle); !errors.Is(err, ErrRunNotFinished) {
		t.Fatalf("early Confirm err=%v, want ErrRunNotFinished", err)
	}
	if err := m.Finish(handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	mustState(t, m, handle, RunStateDraining)
	if err := m.Finish(handle); err != nil {
		t.Fatalf("Finish must be idempotent: %v", err)
	}

	if err := m.Confirm(ctx, handle); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if _, err := m.State(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
	if err := m.Finish(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("Finish after release err=%v, want ErrSessionNotFound", err)
	}
	// released 缓存覆盖重复 Confirm，但 State/Finish 仍按不存在处理。
	if err := m.Confirm(ctx, handle); err != nil {
		t.Fatalf("repeated Confirm err=%v, want nil", err)
	}
	if err := m.Confirm(ctx, ""); !errors.Is(err, ErrEmptySessionID) {
		t.Fatalf("empty Confirm err=%v, want ErrEmptySessionID", err)
	}
}

func TestStartGeneratesUUIDv7(t *testing.T) {
	m := New()
	id, err := m.Start()
	if err != nil {
		t.Fatalf("Start: %v", err)
	}
	parsed, err := uuid.Parse(id)
	if err != nil {
		t.Fatalf("session id %q is not a uuid: %v", id, err)
	}
	if parsed.Version() != 7 {
		t.Fatalf("session id version = %d, want 7", parsed.Version())
	}
}

func TestConfirmWaitsForDrain(t *testing.T) {
	m := New()
	entered := make(chan struct{})
	release := make(chan struct{})

	var s *Session
	handle, err := m.start(sessionOptions{drainFn: func(context.Context) error {
		close(entered)
		<-release
		s.markDrained()
		return nil
	}})
	if err != nil {
		t.Fatalf("start: %v", err)
	}
	s = m.sessions[handle]
	if err := m.Finish(handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}

	done := make(chan error, 1)
	go func() { done <- m.Confirm(context.Background(), handle) }()

	<-entered
	select {
	case err := <-done:
		t.Fatalf("Confirm returned before drain completed: %v", err)
	default:
	}

	close(release)
	if err := <-done; err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if _, err := m.State(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

func TestConcurrentConfirmReleasesOnce(t *testing.T) {
	m := New()
	var releases atomic.Int32
	handle, err := m.start(sessionOptions{releaseFn: func(context.Context) error {
		releases.Add(1)
		return nil
	}})
	if err != nil {
		t.Fatalf("start: %v", err)
	}
	s := m.sessions[handle]
	s.drainFn = func(context.Context) error {
		s.markDrained()
		return nil
	}
	if err := m.Finish(handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}

	const callers = 32
	start := make(chan struct{})
	errs := make(chan error, callers)
	var wg sync.WaitGroup
	for i := 0; i < callers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			<-start
			errs <- m.Confirm(context.Background(), handle)
		}()
	}
	close(start)
	wg.Wait()
	close(errs)
	for err := range errs {
		if err != nil {
			t.Errorf("concurrent Confirm: %v", err)
		}
	}
	if got := releases.Load(); got != 1 {
		t.Fatalf("release hook called %d times, want 1", got)
	}
	// 释放后的重发命中 released 缓存，不得再次执行释放钩子。
	if err := m.Confirm(context.Background(), handle); err != nil {
		t.Fatalf("repeated Confirm: %v", err)
	}
	if got := releases.Load(); got != 1 {
		t.Fatalf("release hook called %d times after retry, want 1", got)
	}
}

func TestConfirmTimeoutThenRetry(t *testing.T) {
	m := New()
	entered := make(chan struct{})
	handle, err := m.start(sessionOptions{drainFn: func(ctx context.Context) error {
		close(entered)
		<-ctx.Done()
		return ctx.Err()
	}})
	if err != nil {
		t.Fatalf("start: %v", err)
	}
	s := m.sessions[handle]
	if err := m.Finish(handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}

	ctx, cancel := context.WithCancel(context.Background())
	done := make(chan error, 1)
	go func() { done <- m.Confirm(ctx, handle) }()
	<-entered
	cancel()
	if err := <-done; !errors.Is(err, context.Canceled) {
		t.Fatalf("Confirm err=%v, want context.Canceled", err)
	}
	mustState(t, m, handle, RunStateDraining)

	// 超时不得污染状态：换掉排空钩子后重试应成功。
	s.drainFn = func(context.Context) error {
		s.markDrained()
		return nil
	}
	if err := m.Confirm(context.Background(), handle); err != nil {
		t.Fatalf("retry Confirm: %v", err)
	}
	if _, err := m.State(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after retry err=%v, want ErrSessionNotFound", err)
	}
}

func TestConfirmFailedSession(t *testing.T) {
	m := New()
	handle, err := m.start(sessionOptions{})
	if err != nil {
		t.Fatalf("start: %v", err)
	}
	s := m.sessions[handle]
	if err := m.Finish(handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	sentinel := errors.New("intake close failed")
	s.fail(sentinel)
	mustState(t, m, handle, RunStateFailed)

	if err := m.Confirm(context.Background(), handle); !errors.Is(err, sentinel) {
		t.Fatalf("Confirm failed session err=%v, want sentinel", err)
	}
	if err := m.Finish(handle); !errors.Is(err, sentinel) {
		t.Fatalf("Finish failed session err=%v, want sentinel", err)
	}
	// FAILED 保留在注册表中，等待后续重试收尾。
	mustState(t, m, handle, RunStateFailed)
}

func TestReleasedCacheEviction(t *testing.T) {
	m := New()
	for i := 0; i < releasedCapacity+4; i++ {
		m.mu.Lock()
		m.released.add(fmt.Sprintf("id-%d", i))
		m.mu.Unlock()
	}
	if m.isReleased("id-0") || m.isReleased("id-3") {
		t.Fatal("oldest released ids must be evicted")
	}
	if !m.isReleased("id-4") || !m.isReleased(fmt.Sprintf("id-%d", releasedCapacity+3)) {
		t.Fatal("recent released ids must stay cached")
	}
	// 重复写入不得占用新的槽位。
	last := fmt.Sprintf("id-%d", releasedCapacity+3)
	m.mu.Lock()
	m.released.add(last)
	m.mu.Unlock()
	if !m.isReleased("id-4") {
		t.Fatal("duplicate add must not evict live entries")
	}
}

func TestRunIsolation(t *testing.T) {
	m := New()
	handleA, err := m.Start()
	if err != nil {
		t.Fatal(err)
	}
	handleB, err := m.Start()
	if err != nil {
		t.Fatal(err)
	}
	if handleA == handleB {
		t.Fatal("two runs must receive distinct handles")
	}
	if err := m.Finish(handleA); err != nil {
		t.Fatal(err)
	}
	mustState(t, m, handleA, RunStateDraining)
	mustState(t, m, handleB, RunStateRunning)
	if err := m.Confirm(context.Background(), handleA); err != nil {
		t.Fatalf("Confirm run A: %v", err)
	}
	mustState(t, m, handleB, RunStateRunning)
}

func TestStartRejectsWhenFull(t *testing.T) {
	m := New()
	handles := make([]string, 0, maxSessions)
	for i := 0; i < maxSessions; i++ {
		handle, err := m.Start()
		if err != nil {
			t.Fatalf("start %d: %v", i, err)
		}
		handles = append(handles, handle)
	}
	if _, err := m.Start(); !errors.Is(err, ErrTooManySessions) {
		t.Fatalf("err = %v, want ErrTooManySessions", err)
	}
	// 释放一个会话后应重新可创建，且不影响其它会话
	if err := m.Finish(handles[0]); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	if err := m.Confirm(context.Background(), handles[0]); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if _, err := m.Start(); err != nil {
		t.Fatalf("Start after release: %v", err)
	}
	mustState(t, m, handles[1], RunStateRunning)
}

func TestConcurrentRuns(t *testing.T) {
	m := New()
	const clients = 4
	var wg sync.WaitGroup
	errs := make(chan error, clients)
	for i := 0; i < clients; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			handle, err := m.Start()
			if err != nil {
				errs <- fmt.Errorf("client %d start: %w", id, err)
				return
			}
			for j := 0; j < 20; j++ {
				if _, err = m.State(handle); err != nil {
					errs <- fmt.Errorf("client %d state: %w", id, err)
					return
				}
			}
			if err = m.Finish(handle); err != nil {
				errs <- fmt.Errorf("client %d finish: %w", id, err)
				return
			}
			if err = m.Confirm(context.Background(), handle); err != nil {
				errs <- fmt.Errorf("client %d confirm: %w", id, err)
			}
		}(i)
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		t.Error(err)
	}
}
