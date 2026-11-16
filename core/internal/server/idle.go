package server

import (
	"context"
	"sync"
	"sync/atomic"
	"time"

	"google.golang.org/grpc/stats"
)

// Timer 抽象空闲定时器句柄，便于测试注入可控定时器。
type Timer interface {
	Stop() bool
}

// AfterFunc 创建一次性定时器，在指定时长后触发回调。
type AfterFunc func(time.Duration, func()) Timer

type stdTimer struct {
	t *time.Timer
}

func (s stdTimer) Stop() bool {
	return s.t.Stop()
}

func defaultAfterFunc(d time.Duration, fn func()) Timer {
	return stdTimer{t: time.AfterFunc(d, fn)}
}

// IdleMonitor 根据活跃连接数维护空闲退出计时。
//
// 当活跃连接数归零且服务已启动空闲监控时开始计时；新连接到达立即取消计时。
// 超时后最多触发一次 onIdle 回调，且在持锁外执行以避免死锁。
type IdleMonitor struct {
	timeout   time.Duration
	onIdle    func()
	afterFunc AfterFunc

	mu      sync.Mutex
	active  int64
	armed   bool
	stopped bool
	gen     uint64
	timer   Timer
}

// IdleOption 配置 IdleMonitor 的可选行为。
type IdleOption func(*IdleMonitor)

// WithAfterFunc 注入自定义定时器工厂，供单测替代真实时间等待。
func WithAfterFunc(fn AfterFunc) IdleOption {
	return func(m *IdleMonitor) {
		if fn != nil {
			m.afterFunc = fn
		}
	}
}

// NewIdleMonitor 创建空闲监控器；timeout <= 0 时不会启动空闲计时。
func NewIdleMonitor(timeout time.Duration, onIdle func(), opts ...IdleOption) *IdleMonitor {
	m := &IdleMonitor{
		timeout:   timeout,
		onIdle:    onIdle,
		afterFunc: defaultAfterFunc,
	}
	for _, opt := range opts {
		if opt != nil {
			opt(m)
		}
	}
	return m
}

// Start 开启空闲监控；若此时没有活跃连接，立即进入空闲倒计时。
func (m *IdleMonitor) Start() {
	m.mu.Lock()
	defer m.mu.Unlock()
	if m.stopped || m.timeout <= 0 {
		return
	}
	m.armed = true
	if m.active == 0 && m.timer == nil {
		m.armLocked()
	}
}

// Stop 永久停止空闲监控并取消当前定时器，可重复调用。
func (m *IdleMonitor) Stop() {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.stopped = true
	m.gen++
	if m.timer != nil {
		m.timer.Stop()
		m.timer = nil
	}
}

// OnConnBegin 记录一条新连接进入，并取消尚未触发的空闲计时。
func (m *IdleMonitor) OnConnBegin() {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.active++
	m.gen++
	if m.timer != nil {
		m.timer.Stop()
		m.timer = nil
	}
}

// OnConnEnd 记录一条连接断开；连接数归零时重新开始空闲倒计时。
func (m *IdleMonitor) OnConnEnd() {
	m.mu.Lock()
	defer m.mu.Unlock()
	if m.active > 0 {
		m.active--
	}
	m.gen++
	if m.armed && !m.stopped && m.active == 0 && m.timer == nil {
		m.armLocked()
	}
}

// Active 返回当前记录的活跃连接数。
func (m *IdleMonitor) Active() int64 {
	m.mu.Lock()
	defer m.mu.Unlock()
	return m.active
}

func (m *IdleMonitor) armLocked() {
	myGen := m.gen
	m.timer = m.afterFunc(m.timeout, func() {
		m.fire(myGen)
	})
}

func (m *IdleMonitor) fire(myGen uint64) {
	m.mu.Lock()
	if m.stopped || !m.armed || m.active != 0 || m.gen != myGen {
		m.mu.Unlock()
		return
	}
	m.stopped = true
	m.timer = nil
	cb := m.onIdle
	m.mu.Unlock()

	if cb != nil {
		cb()
	}
}

type connCtxKey struct{}

type connState struct {
	closed atomic.Bool
}

// ConnTracker 实现 gRPC stats.Handler，用于统计传输层活跃连接并驱动 IdleMonitor。
type ConnTracker struct {
	active atomic.Int64

	mu      sync.RWMutex
	monitor *IdleMonitor
}

// NewConnTracker 创建连接追踪器。
func NewConnTracker() *ConnTracker {
	return &ConnTracker{}
}

// SetMonitor 绑定空闲监控器；应在 gRPC Serve 开始前设置。
func (t *ConnTracker) SetMonitor(m *IdleMonitor) {
	t.mu.Lock()
	t.monitor = m
	t.mu.Unlock()
}

// Active 返回当前活跃传输连接数。
func (t *ConnTracker) Active() int64 {
	return t.active.Load()
}

// TagRPC 满足 stats.Handler 接口，连接追踪无需附加 RPC 上下文。
func (t *ConnTracker) TagRPC(ctx context.Context, _ *stats.RPCTagInfo) context.Context {
	return ctx
}

// HandleRPC 满足 stats.Handler 接口，连接追踪只关心连接生命周期事件。
func (t *ConnTracker) HandleRPC(context.Context, stats.RPCStats) {}

// TagConn 为每条传输连接绑定幂等关闭标记。
func (t *ConnTracker) TagConn(ctx context.Context, _ *stats.ConnTagInfo) context.Context {
	return context.WithValue(ctx, connCtxKey{}, &connState{})
}

// HandleConn 统计连接建立与断开事件，并同步通知绑定的 IdleMonitor。
func (t *ConnTracker) HandleConn(ctx context.Context, st stats.ConnStats) {
	switch st.(type) {
	case *stats.ConnBegin:
		t.active.Add(1)
		if m := t.getMonitor(); m != nil {
			m.OnConnBegin()
		}
	case *stats.ConnEnd:
		if state, ok := ctx.Value(connCtxKey{}).(*connState); ok {
			if !state.closed.CompareAndSwap(false, true) {
				return
			}
		}
		t.active.Add(-1)
		if m := t.getMonitor(); m != nil {
			m.OnConnEnd()
		}
	}
}

func (t *ConnTracker) getMonitor() *IdleMonitor {
	t.mu.RLock()
	defer t.mu.RUnlock()
	return t.monitor
}
