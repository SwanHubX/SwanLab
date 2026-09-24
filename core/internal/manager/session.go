// Package manager 负责维护 swanlab-core 的 Run 会话与业务编排逻辑。
//
// 该包与具体传输协议无关（不依赖 gRPC/Protobuf），专职管理 Run 级状态机与会话路由。
package manager

import "sync"

// RunState 表示单个实验会话的内部生命周期状态。
type RunState int32

const (
	// RunStateRunning 会话正在运行中。
	RunStateRunning RunState = iota
	// RunStateFinished 会话已标记完成。
	RunStateFinished
)

// Session 保存单个 run_handle 对应的会话状态与生命周期数据。
type Session struct {
	mu    sync.Mutex
	state RunState
}

func newSession() *Session {
	return &Session{state: RunStateRunning}
}

// finish 幂等标记会话已结束。
func (s *Session) finish() {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.state = RunStateFinished
}

// stateSnapshot 返回当前会话状态的并发安全快照。
func (s *Session) stateSnapshot() RunState {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.state
}

// finished 检查会话是否已完成。
func (s *Session) finished() bool {
	return s.stateSnapshot() == RunStateFinished
}
