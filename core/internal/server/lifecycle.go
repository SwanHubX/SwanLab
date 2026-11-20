package server

import "sync"

// ServiceState 表示 gRPC 服务的进程级生命周期状态。
type ServiceState int32

const (
	// StateNotReady 服务已监听但尚未完成握手（SpinupService），拒绝业务 RPC。
	StateNotReady ServiceState = iota
	// StateReady 服务已就绪，可正常接受和处理业务请求。
	StateReady
	// StateStopping 服务正在执行退出流程。
	StateStopping
	// StateClosed 服务已完全关闭。
	StateClosed
)

// Lifecycle 维护并发安全的服务生命周期状态机。
//
// 状态流转：
// - StateNotReady -> StateReady：通过 Spinup 触发，Ready 状态下重复调用幂等成功；
// - StateNotReady / StateReady -> StateStopping：服务开始关闭（Shutdown）时进入；
// - StateStopping -> StateClosed：服务收尾完成后进入。
type Lifecycle struct {
	mu    sync.RWMutex
	state ServiceState
}

func newLifecycle() *Lifecycle {
	return &Lifecycle{state: StateNotReady}
}

// Get 返回当前的服务生命周期状态。
func (l *Lifecycle) Get() ServiceState {
	l.mu.RLock()
	defer l.mu.RUnlock()
	return l.state
}

// Ready 检查服务当前是否处于就绪（StateReady）状态。
func (l *Lifecycle) Ready() bool {
	return l.Get() == StateReady
}

// Spinup 将服务状态置为 StateReady。
// 在 StateReady 状态下重复调用幂等返回 true；若服务已处于退出或关闭中，则返回 false。
func (l *Lifecycle) Spinup() bool {
	l.mu.Lock()
	defer l.mu.Unlock()
	switch l.state {
	case StateNotReady:
		l.state = StateReady
		return true
	case StateReady:
		return true
	default:
		return false
	}
}

// BeginStopping 将服务置为 StateStopping 退出状态，多次调用幂等。
func (l *Lifecycle) BeginStopping() {
	l.mu.Lock()
	defer l.mu.Unlock()
	if l.state == StateNotReady || l.state == StateReady {
		l.state = StateStopping
	}
}

// Close 将服务置为 StateClosed 状态，标识服务已完全关闭。
func (l *Lifecycle) Close() {
	l.mu.Lock()
	defer l.mu.Unlock()
	l.state = StateClosed
}
