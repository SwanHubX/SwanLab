package health

import (
	"context"
	"sync"

	"github.com/swanhubx/swanlab/core/internal/cron"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// Site 聚合一个站点的控制面与数据面：任一平面连续不可达达到阈值即降级。
// 写入侧降级通过 Report(ReasonWrite, ...) 注入，与探针状态取或。
type Site struct {
	key     string
	control *prober
	data    *prober

	updates chan struct{}

	mu           sync.Mutex
	closed       bool
	probeControl bool
	probeData    bool
	write        bool
	state        State
	reason       Reason

	stopOnce  sync.Once
	scheduler *cron.Scheduler
}

func newSite(opts SiteOptions, clock cron.Clock) *Site {
	interval := opts.Interval
	if interval <= 0 {
		interval = DefaultInterval
	}
	timeout := opts.Timeout
	if timeout <= 0 {
		timeout = DefaultTimeout
	}
	degradeAfter := opts.DegradeThreshold
	if degradeAfter <= 0 {
		degradeAfter = DefaultDegradeThreshold
	}
	recoverAfter := opts.RecoverThreshold
	if recoverAfter <= 0 {
		recoverAfter = DefaultRecoverThreshold
	}

	s := &Site{key: opts.Key, updates: make(chan struct{}, 1)}
	probeOpts := proberOptions{
		timeout:      timeout,
		degradeAfter: degradeAfter,
		recoverAfter: recoverAfter,
	}
	s.control = newProber(opts.Control, probeOpts, func(degraded bool) {
		s.setReason(ReasonProbeControl, degraded)
	})
	s.data = newProber(opts.Data, probeOpts, func(degraded bool) {
		s.setReason(ReasonProbeData, degraded)
	})
	// 启动预检在任务协程启动前写入初值，避免首个探测周期内的错误乐观状态。
	if opts.InitialControlDegraded {
		s.control.seedDegraded()
		s.setReason(ReasonProbeControl, true)
	}
	if opts.InitialDataDegraded {
		s.data.seedDegraded()
		s.setReason(ReasonProbeData, true)
	}

	s.scheduler = cron.New(clock)
	s.scheduler.Add(cron.Job{Name: "control", Interval: interval, Run: s.control.tick})
	s.scheduler.Add(cron.Job{Name: "data", Interval: interval, Run: s.data.tick})
	s.scheduler.Start(context.Background())
	return s
}

// Key 返回站点标识（api_host + api_key 摘要），同一 Key 在 Registry 内共享实例。
func (s *Site) Key() string { return s.key }

// Updates 在聚合状态或降级原因变化时收到通知，容量 1，多次变化会合并。
// 消费方应配合 Snapshot 读取详情；Site 释放后该通道关闭。
func (s *Site) Updates() <-chan struct{} { return s.updates }

// MergeInitial 把后续会话的启动预检结论并入共享站点；只降级、不升级。
// 探针随后仍按连续成功阈值恢复该平面。
func (s *Site) MergeInitial(controlDegraded, dataDegraded bool) {
	if controlDegraded {
		s.control.degradeNow()
	}
	if dataDegraded {
		s.data.degradeNow()
	}
}

// Report 注入或解除写入侧降级信号，是否解除由上层按“上传成功一次”判定。
// 探针原因由 prober 自己维护；外部注入探针原因会与阈值状态脱节，因此直接 panic。
func (s *Site) Report(reason Reason, active bool) {
	if reason != ReasonWrite {
		panic("health: only ReasonWrite can be reported externally")
	}
	s.setReason(reason, active)
}

// Snapshot 返回当前聚合状态与两个平面的最近探测结果。
func (s *Site) Snapshot() Snapshot {
	s.mu.Lock()
	state, reason := s.state, s.reason
	s.mu.Unlock()
	return Snapshot{
		State:            state,
		Reason:           reason,
		ControlReachable: s.control.reachable.Load(),
		DataReachable:    s.data.reachable.Load(),
	}
}

func (s *Site) setReason(reason Reason, active bool) {
	s.mu.Lock()
	if s.closed {
		s.mu.Unlock()
		return
	}
	switch reason {
	case ReasonProbeControl:
		s.probeControl = active
	case ReasonProbeData:
		s.probeData = active
	case ReasonWrite:
		s.write = active
	default:
		s.mu.Unlock()
		return
	}
	state, next := s.aggregateLocked()
	if state == s.state && next == s.reason {
		s.mu.Unlock()
		return
	}
	s.state, s.reason = state, next
	select {
	case s.updates <- struct{}{}:
	default:
	}
	s.mu.Unlock()

	// 状态翻转只记录一次，周期性探测结果不刷屏。
	if state == StateDegraded {
		console.Warningf("site health degraded: reason=%s", next)
	} else {
		console.Infof("site health recovered: state=%s", state)
	}
}

// aggregateLocked 按控制面 → 数据面 → 写入的优先级给出唯一原因。
// 控制面在数据面关键路径上（列定义先于标量 / 媒体），因此优先展示。
func (s *Site) aggregateLocked() (State, Reason) {
	switch {
	case s.probeControl:
		return StateDegraded, ReasonProbeControl
	case s.probeData:
		return StateDegraded, ReasonProbeData
	case s.write:
		return StateDegraded, ReasonWrite
	default:
		return StateNormal, ""
	}
}

// stop 停止探测任务并关闭通知通道；可重复调用。
func (s *Site) stop() {
	s.stopOnce.Do(func() {
		s.scheduler.Stop()
		s.mu.Lock()
		s.closed = true
		close(s.updates)
		s.mu.Unlock()
	})
}
