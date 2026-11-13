package manager

import (
	"context"
	"errors"
	"fmt"
	"strings"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"github.com/google/uuid"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
)

// drainBare 直接驱动裸会话进入 DRAINING，用于不装配收尾资源的状态机测试。
// finishCommitted 置位以模拟 FinishRecord 已提交，避免触发“Shutdown 仅排空”的保护。
func drainBare(m *Manager, handle string) error {
	s, err := m.session(handle)
	if err != nil {
		return err
	}
	s.finishCommitted = true
	return s.beginDrain()
}

// confirm 只关心 Confirm 的错误，供不校验响应文案的状态机测试使用。
func confirm(ctx context.Context, m *Manager, handle string) error {
	_, err := m.Confirm(ctx, handle)
	return err
}

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
	m := newTestManager()
	ctx := context.Background()
	handle := installBare(t, m, sessionOptions{})
	mustState(t, m, handle, RunStateRunning)

	if err := confirm(ctx, m, handle); !errors.Is(err, ErrRunNotFinished) {
		t.Fatalf("early Confirm err=%v, want ErrRunNotFinished", err)
	}
	if err := drainBare(m, handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	mustState(t, m, handle, RunStateDraining)
	if err := drainBare(m, handle); err != nil {
		t.Fatalf("Finish must be idempotent: %v", err)
	}

	if err := confirm(ctx, m, handle); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if _, err := m.State(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
	if err := drainBare(m, handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("Finish after release err=%v, want ErrSessionNotFound", err)
	}
	// released 缓存覆盖重复 Confirm，但 State/Finish 仍按不存在处理。
	if err := confirm(ctx, m, handle); err != nil {
		t.Fatalf("repeated Confirm err=%v, want nil", err)
	}
	if err := confirm(ctx, m, ""); !errors.Is(err, ErrEmptySessionID) {
		t.Fatalf("empty Confirm err=%v, want ErrEmptySessionID", err)
	}
}

func TestInstallGeneratesUUIDv7(t *testing.T) {
	m := newTestManager()
	handle := installBare(t, m, sessionOptions{})
	parsed, err := uuid.Parse(handle)
	if err != nil {
		t.Fatalf("session id %q is not a uuid: %v", handle, err)
	}
	if parsed.Version() != 7 {
		t.Fatalf("session id version = %d, want 7", parsed.Version())
	}
}

func TestConfirmWaitsForDrain(t *testing.T) {
	m := newTestManager()
	entered := make(chan struct{})
	release := make(chan struct{})
	handle := installBare(t, m, sessionOptions{})
	s := m.sessions[handle]
	s.drainFn = func(context.Context) error {
		close(entered)
		<-release
		s.markDrained()
		return nil
	}
	if err := drainBare(m, handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}

	done := make(chan error, 1)
	go func() { done <- confirm(context.Background(), m, handle) }()

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
	m := newTestManager()
	var releases atomic.Int32
	handle := installBare(t, m, sessionOptions{releaseFn: func(context.Context) error {
		releases.Add(1)
		return nil
	}})
	s := m.sessions[handle]
	s.drainFn = func(context.Context) error {
		s.markDrained()
		return nil
	}
	if err := drainBare(m, handle); err != nil {
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
			errs <- confirm(context.Background(), m, handle)
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
	if err := confirm(context.Background(), m, handle); err != nil {
		t.Fatalf("repeated Confirm: %v", err)
	}
	if got := releases.Load(); got != 1 {
		t.Fatalf("release hook called %d times after retry, want 1", got)
	}
}

// 收尾失败必须能唤醒正在等排空的 Confirm：否则 drainFn 永不返回，
// Confirm 只能等到自己的 deadline，把可重试的失败变成超时。
func TestConfirmWakesOnFinishFailure(t *testing.T) {
	m := newTestManager()
	entered := make(chan struct{})
	handle := installBare(t, m, sessionOptions{})
	s := m.sessions[handle]
	s.drainFn = func(ctx context.Context) error {
		close(entered)
		<-ctx.Done()
		return ctx.Err()
	}
	if err := drainBare(m, handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}

	done := make(chan error, 1)
	go func() { done <- confirm(context.Background(), m, handle) }()
	<-entered
	s.fail(errors.New("boom"))

	select {
	case err := <-done:
		if !errors.Is(err, ErrRunFailed) {
			t.Fatalf("Confirm err=%v, want ErrRunFailed", err)
		}
		if errors.Is(err, context.Canceled) {
			t.Fatalf("Confirm err=%v, must not be reported as cancellation", err)
		}
	case <-time.After(2 * time.Second):
		t.Fatal("Confirm did not wake after the finish failed")
	}
}

func TestConfirmTimeoutThenRetry(t *testing.T) {
	m := newTestManager()
	entered := make(chan struct{})
	handle := installBare(t, m, sessionOptions{})
	s := m.sessions[handle]
	s.drainFn = func(ctx context.Context) error {
		close(entered)
		<-ctx.Done()
		return ctx.Err()
	}
	if err := drainBare(m, handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}

	ctx, cancel := context.WithCancel(context.Background())
	done := make(chan error, 1)
	go func() { done <- confirm(ctx, m, handle) }()
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
	if err := confirm(context.Background(), m, handle); err != nil {
		t.Fatalf("retry Confirm: %v", err)
	}
	if _, err := m.State(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after retry err=%v, want ErrSessionNotFound", err)
	}
}

func TestConfirmFailedSession(t *testing.T) {
	m := newTestManager()
	handle := installBare(t, m, sessionOptions{})
	s := m.sessions[handle]
	if err := drainBare(m, handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	sentinel := errors.New("intake close failed")
	s.fail(sentinel)
	mustState(t, m, handle, RunStateFailed)

	if err := confirm(context.Background(), m, handle); !errors.Is(err, sentinel) {
		t.Fatalf("Confirm failed session err=%v, want sentinel", err)
	}
	if err := drainBare(m, handle); !errors.Is(err, sentinel) {
		t.Fatalf("Finish failed session err=%v, want sentinel", err)
	}
	// FAILED 保留在注册表中，等待后续重试收尾。
	mustState(t, m, handle, RunStateFailed)
}

// stop_experiment 失败必须保留会话：第一次返回 release-failed 响应，
// 第二次 Confirm 重试成功后才移除，且释放钩子恰好被调用两次。
func TestConfirmReleaseFailureKeepsSessionForRetry(t *testing.T) {
	m := newTestManager()
	var attempts atomic.Int32
	handle := installBare(t, m, sessionOptions{releaseFn: func(context.Context) error {
		if attempts.Add(1) == 1 {
			return errors.New("stop_experiment unavailable")
		}
		return nil
	}})
	s := m.sessions[handle]
	s.drainFn = func(context.Context) error {
		s.markDrained()
		return nil
	}
	if err := drainBare(m, handle); err != nil {
		t.Fatalf("Finish: %v", err)
	}

	result, err := m.Confirm(context.Background(), handle)
	if err != nil {
		t.Fatalf("first Confirm err=%v", err)
	}
	if !result.ReleaseFailed {
		t.Fatalf("first Confirm result=%+v, want release failure", result)
	}
	if !strings.Contains(result.Message, "saved locally") {
		t.Fatalf("first Confirm message=%q, want local-copy hint", result.Message)
	}
	mustState(t, m, handle, RunStateDrained)

	result, err = m.Confirm(context.Background(), handle)
	if err != nil {
		t.Fatalf("retry Confirm err=%v", err)
	}
	if result.ReleaseFailed {
		t.Fatalf("retry Confirm result=%+v, want success", result)
	}
	if got := attempts.Load(); got != 2 {
		t.Fatalf("release hook called %d times, want 2", got)
	}
	if _, err = m.State(handle); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after retry err=%v, want ErrSessionNotFound", err)
	}
}

func TestConfirmReleaseFailureClassification(t *testing.T) {
	cases := []struct {
		name string
		err  error
		want bool
	}{
		{"network unavailable", &client.APIError{Retryable: true}, true},
		{"deadline exceeded", context.DeadlineExceeded, true},
		{"deterministic 4xx", &client.APIError{StatusCode: 400}, false},
		{"plain error", errors.New("boom"), false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			m := newTestManager()
			handle := installBare(t, m, sessionOptions{releaseFn: func(context.Context) error {
				return tc.err
			}})
			s := m.sessions[handle]
			s.drainFn = func(context.Context) error {
				s.markDrained()
				return nil
			}
			if err := drainBare(m, handle); err != nil {
				t.Fatalf("Finish: %v", err)
			}

			result, err := m.Confirm(context.Background(), handle)
			if err != nil {
				t.Fatalf("Confirm err=%v", err)
			}
			if !result.ReleaseFailed {
				t.Fatalf("result=%+v, want release failure", result)
			}
			if result.Retryable != tc.want {
				t.Fatalf("Retryable=%v, want %v", result.Retryable, tc.want)
			}
			mustState(t, m, handle, RunStateDrained)
		})
	}
}

func TestReleasedCacheEviction(t *testing.T) {
	m := newTestManager()
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
	m := newTestManager()
	handleA := installBare(t, m, sessionOptions{})
	handleB := installBare(t, m, sessionOptions{})
	if handleA == handleB {
		t.Fatal("two runs must receive distinct handles")
	}
	if err := drainBare(m, handleA); err != nil {
		t.Fatal(err)
	}
	mustState(t, m, handleA, RunStateDraining)
	mustState(t, m, handleB, RunStateRunning)
	if err := confirm(context.Background(), m, handleA); err != nil {
		t.Fatalf("Confirm run A: %v", err)
	}
	mustState(t, m, handleB, RunStateRunning)
}

func TestStartRejectsWhenFull(t *testing.T) {
	m := newTestManager()
	handles := make([]string, 0, maxSessions)
	for i := 0; i < maxSessions; i++ {
		handles = append(handles, installBare(t, m, sessionOptions{}))
	}
	if _, err := m.install(newSession(sessionOptions{})); !errors.Is(err, ErrTooManySessions) {
		t.Fatalf("err = %v, want ErrTooManySessions", err)
	}
	// 释放一个会话后应重新可创建，且不影响其它会话
	if err := drainBare(m, handles[0]); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	if err := confirm(context.Background(), m, handles[0]); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if _, err := m.install(newSession(sessionOptions{})); err != nil {
		t.Fatalf("install after release: %v", err)
	}
	mustState(t, m, handles[1], RunStateRunning)
}

func TestConcurrentRuns(t *testing.T) {
	m := newTestManager()
	const clients = 4
	var wg sync.WaitGroup
	errs := make(chan error, clients)
	for i := 0; i < clients; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			handle, err := m.install(newSession(sessionOptions{}))
			if err != nil {
				errs <- fmt.Errorf("client %d install: %w", id, err)
				return
			}
			for j := 0; j < 20; j++ {
				if _, err = m.State(handle); err != nil {
					errs <- fmt.Errorf("client %d state: %w", id, err)
					return
				}
			}
			if err = drainBare(m, handle); err != nil {
				errs <- fmt.Errorf("client %d finish: %w", id, err)
				return
			}
			if err = confirm(context.Background(), m, handle); err != nil {
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

// Shutdown 面向进程退出：排空会话、刷盘且不写 FinishRecord，也不上报实验终态。
func TestManagerShutdownDrainsActiveSessions(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("shutdown/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	var stops int
	h.session.releaseFn = func(context.Context) error {
		stops++
		return nil
	}
	if err := h.mgr.Shutdown(ctx); err != nil {
		t.Fatalf("Shutdown: %v", err)
	}
	if stops != 0 {
		t.Fatalf("Shutdown must not call stop_experiment, got %d call(s)", stops)
	}
	if state := h.session.stateSnapshot(); state != RunStateDrained {
		t.Fatalf("state = %s, want drained", state)
	}

	records := readStore(t, h.runFile)
	if len(records) == 0 {
		t.Fatal("shutdown must flush accepted records to the local store")
	}
	for _, item := range records {
		if item.GetNum() == record.NumFinish {
			t.Fatal("shutdown must not synthesize a FinishRecord")
		}
	}
	if got := len(h.uploads.snapshot()); got == 0 {
		t.Fatal("shutdown must drain accepted records through the uploader")
	}
}

// Shutdown 只做上传排空、不写 FinishRecord：Confirm 不能把它当作正常结束释放。
func TestConfirmRejectsSessionDrainedByShutdown(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("shutdown/confirm", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	var stops int
	h.session.releaseFn = func(context.Context) error {
		stops++
		return nil
	}

	if err := h.mgr.Shutdown(ctx); err != nil {
		t.Fatalf("Shutdown: %v", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateDrained {
		t.Fatalf("state after shutdown = %s, want drained", state)
	}

	if _, err := h.mgr.Confirm(ctx, h.id); !errors.Is(err, ErrRunFailed) {
		t.Fatalf("Confirm err=%v, want ErrRunFailed", err)
	}
	if stops != 0 {
		t.Fatalf("stop_experiment called %d times, want 0", stops)
	}
	if _, err := h.mgr.State(h.id); err != nil {
		t.Fatalf("session must stay registered after rejected Confirm: %v", err)
	}
	if h.mgr.isReleased(h.id) {
		t.Fatal("rejected session must not enter the released cache")
	}
}
