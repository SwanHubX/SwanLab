// Package manager 负责维护 swanlab-core 的 Run 会话与业务编排逻辑。
//
// 该包与具体传输协议无关（不依赖 gRPC 服务契约），专职管理 Run 级状态机与会话路由。
package manager

import (
	"context"
	"sync"
	"sync/atomic"
)

// RunState 表示单个实验会话的内部生命周期状态。
type RunState int32

const (
	// RunStateRunning 会话正在运行中，接受公共提交。
	RunStateRunning RunState = iota
	// RunStateDraining 会话已停止接受公共提交，正在等待上传与落盘排空。
	RunStateDraining
	// RunStateDrained 会话已排空，可以安全释放。
	RunStateDrained
	// RunStateFailed 收尾序列失败，会话不能按成功路径释放。
	RunStateFailed
)

// String 返回状态的可读名称，便于日志与测试失败信息定位。
func (s RunState) String() string {
	switch s {
	case RunStateRunning:
		return "running"
	case RunStateDraining:
		return "draining"
	case RunStateDrained:
		return "drained"
	case RunStateFailed:
		return "failed"
	default:
		return "unknown"
	}
}

// sessionOptions 注入会话的排空与释放钩子。
// C1 骨架允许为空；C3 会把 drainFn 接到 uploader.Wait，C6 把 releaseFn 接到 stop_experiment。
type sessionOptions struct {
	drainFn   func(context.Context) error
	releaseFn func(context.Context) error
}

// Session 保存单个 session_id 对应的会话状态与生命周期数据。
type Session struct {
	// mu 保护状态切换与失败信息，与公共提交共用同一把锁。
	mu sync.Mutex
	// confirmMu 只串行化「释放」阶段，不覆盖排空等待。
	confirmMu sync.Mutex

	state atomic.Int32

	// drained 在 DRAINING → DRAINED 时关闭；failed 在进入 FAILED 时关闭。
	drained chan struct{}
	failed  chan struct{}

	drainFn   func(context.Context) error
	releaseFn func(context.Context) error

	finishErr error
}

func newSession(opts sessionOptions) *Session {
	s := &Session{
		drained:   make(chan struct{}),
		failed:    make(chan struct{}),
		drainFn:   opts.drainFn,
		releaseFn: opts.releaseFn,
	}
	s.state.Store(int32(RunStateRunning))
	return s
}

// beginDrain 把会话从 RUNNING 切到 DRAINING；重复调用幂等。
// 已 FAILED 时返回原始失败，调用方不得把它当作成功。
func (s *Session) beginDrain() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	switch RunState(s.state.Load()) {
	case RunStateRunning:
		s.state.Store(int32(RunStateDraining))
		return nil
	case RunStateDraining, RunStateDrained:
		return nil
	case RunStateFailed:
		return s.failureLocked()
	default:
		return ErrRunFailed
	}
}

// markDrained 由排空完成回调调用；只允许 DRAINING → DRAINED。
func (s *Session) markDrained() bool {
	if !s.state.CompareAndSwap(int32(RunStateDraining), int32(RunStateDrained)) {
		return false
	}
	close(s.drained)
	return true
}

// fail 记录收尾失败并把会话置为 FAILED；仅对 DRAINING 生效。
func (s *Session) fail(err error) {
	if err == nil {
		err = ErrRunFailed
	}
	s.mu.Lock()
	defer s.mu.Unlock()
	if RunState(s.state.Load()) != RunStateDraining {
		return
	}
	s.finishErr = err
	s.state.Store(int32(RunStateFailed))
	close(s.failed)
}

// waitDrained 在锁外等待会话排空，返回 nil 表示可以进入释放阶段。
func (s *Session) waitDrained(ctx context.Context) error {
	switch RunState(s.state.Load()) {
	case RunStateRunning:
		return ErrRunNotFinished
	case RunStateDrained:
		return nil
	case RunStateFailed:
		return s.failure()
	}

	// C1 骨架尚未装配异步排空器；没有待排空数据时直接收敛，
	// C3 注入 uploader.Wait 后走下面的等待分支。
	if s.drainFn == nil {
		s.markDrained()
		return nil
	}
	if err := s.drainFn(ctx); err != nil {
		select {
		case <-s.failed:
			return s.failure()
		default:
			return err
		}
	}
	select {
	case <-s.drained:
		return nil
	case <-s.failed:
		return s.failure()
	case <-ctx.Done():
		return ctx.Err()
	}
}

// stateSnapshot 返回当前会话状态的并发安全快照。
func (s *Session) stateSnapshot() RunState {
	return RunState(s.state.Load())
}

// failure 返回收尾失败信息。
func (s *Session) failure() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.failureLocked()
}

func (s *Session) failureLocked() error {
	if s.finishErr != nil {
		return s.finishErr
	}
	return ErrRunFailed
}
