package server

import "sync"

// ServiceState 是服务级生命周期状态，只存在于 Go 内存，不上 proto。
// 线上表达是 gRPC 状态码：READY 之前的 run RPC 返回 FAILED_PRECONDITION，
// 关闭仲裁由 Controller 汇入统一关闭路径。
type ServiceState int32

const (
	// StateNotReady 服务已监听但未完成 Spinup，拒绝 run 级 RPC。
	StateNotReady ServiceState = iota
	// StateReady SpinupService 成功，接受 run 级 RPC。
	StateReady
	// StateStopping 关闭序列进行中（Teardown/信号/父死/serve 错误）。
	StateStopping
	// StateClosed 关闭序列完成。
	StateClosed
)

// Lifecycle 是服务状态机的并发安全实现。
//
// 转换：NOT_READY → READY（SpinupService，READY 下幂等）；
// NOT_READY/READY → STOPPING（Controller.Shutdown，幂等）；STOPPING → CLOSED（关闭完成）。
type Lifecycle struct {
	mu    sync.RWMutex
	state ServiceState
}

func newLifecycle() *Lifecycle {
	return &Lifecycle{state: StateNotReady}
}

// Get 返回当前服务状态。
func (l *Lifecycle) Get() ServiceState {
	l.mu.RLock()
	defer l.mu.RUnlock()
	return l.state
}

// Ready 报告服务是否处于 READY。
func (l *Lifecycle) Ready() bool {
	return l.Get() == StateReady
}

// Spinup 幂等执行 NOT_READY → READY；READY 下重复调用成功；
// 已进入 STOPPING/CLOSED 时返回 false，由调用方映射为 FAILED_PRECONDITION。
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

// BeginStopping 幂等进入 STOPPING；已在 STOPPING/CLOSED 时为 no-op。
func (l *Lifecycle) BeginStopping() {
	l.mu.Lock()
	defer l.mu.Unlock()
	if l.state == StateNotReady || l.state == StateReady {
		l.state = StateStopping
	}
}

// Close 在关闭序列完成后进入 CLOSED；幂等。
func (l *Lifecycle) Close() {
	l.mu.Lock()
	defer l.mu.Unlock()
	l.state = StateClosed
}
