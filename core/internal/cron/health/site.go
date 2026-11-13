package health

import (
	"context"
	"sync"

	"github.com/swanhubx/swanlab/core/internal/cron"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// Site 聚合一个站点的 swanlab-server 与 swanlab-house：任一服务连续不可达达到阈值即降级。
// 写入侧降级通过 Report(ReasonWrite, ...) 注入，与探针状态取或。
type Site struct {
	key    string
	server *prober
	house  *prober

	mu             sync.Mutex
	closed         bool
	serverDegraded bool
	houseDegraded  bool
	write          bool
	state          State
	reason         Reason

	nextSubscription uint64
	subscriptions    map[uint64]chan struct{}

	stopOnce  sync.Once
	scheduler *cron.Scheduler
}

func newSite(opts SiteOptions, clock cron.Clock) *Site {
	interval := positiveOrDefault(opts.Interval, DefaultInterval)
	timeout := positiveOrDefault(opts.Timeout, DefaultTimeout)
	degradeAfter := positiveOrDefault(opts.DegradeThreshold, DefaultDegradeThreshold)
	recoverAfter := positiveOrDefault(opts.RecoverThreshold, DefaultRecoverThreshold)

	s := &Site{key: opts.Key, subscriptions: make(map[uint64]chan struct{})}
	probeOpts := proberOptions{
		timeout:      timeout,
		degradeAfter: degradeAfter,
		recoverAfter: recoverAfter,
	}
	s.server = newProber(opts.Server, probeOpts, func(degraded bool) {
		s.setReason(ReasonProbeServer, degraded)
	})
	s.house = newProber(opts.House, probeOpts, func(degraded bool) {
		s.setReason(ReasonProbeHouse, degraded)
	})
	// 启动预检在任务协程启动前写入初值，避免首个探测周期内的错误乐观状态。
	if opts.InitialHouseDegraded {
		s.house.degradeNow()
	}

	s.scheduler = cron.New(clock)
	s.scheduler.Add(cron.Job{Interval: interval, Run: s.server.tick})
	s.scheduler.Add(cron.Job{Interval: interval, Run: s.house.tick})
	s.scheduler.Start(context.Background())
	return s
}

// Subscribe 返回独立的状态通知 channel 与幂等注销函数。
// 状态变化会非阻塞通知所有订阅者；消费方应配合 Snapshot 读取最新状态。
func (s *Site) Subscribe() (<-chan struct{}, func()) {
	s.mu.Lock()
	defer s.mu.Unlock()
	updates := make(chan struct{}, 1)
	if s.closed {
		close(updates)
		return updates, func() {}
	}
	id := s.nextSubscription
	s.nextSubscription++
	s.subscriptions[id] = updates
	unsubscribe := func() {
		s.mu.Lock()
		defer s.mu.Unlock()
		if current, ok := s.subscriptions[id]; ok {
			delete(s.subscriptions, id)
			close(current)
		}
	}
	return updates, unsubscribe
}

// MergeInitial 把后续会话的启动预检结论并入共享站点；只降级、不升级。
// 探针随后仍按连续成功阈值恢复该服务。
func (s *Site) MergeInitial(houseDegraded bool) {
	if houseDegraded {
		s.house.degradeNow()
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

// Snapshot 返回当前聚合状态与降级来源。
func (s *Site) Snapshot() Snapshot {
	s.mu.Lock()
	defer s.mu.Unlock()
	return Snapshot{State: s.state, Reason: s.reason}
}

func (s *Site) setReason(reason Reason, active bool) {
	s.mu.Lock()
	if s.closed {
		s.mu.Unlock()
		return
	}
	switch reason {
	case ReasonProbeServer:
		s.serverDegraded = active
	case ReasonProbeHouse:
		s.houseDegraded = active
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
	for _, updates := range s.subscriptions {
		select {
		case updates <- struct{}{}:
		default:
		}
	}
	s.mu.Unlock()

	// 状态翻转只记录一次，周期性探测结果不刷屏。
	if state == StateDegraded {
		console.Warningf("site health degraded: reason=%s", next)
	} else {
		console.Infof("site health recovered: state=%s", state)
	}
}

// aggregateLocked 按 server → house → 写入的优先级给出唯一原因。
// server 在 house 关键路径上（列定义先于标量 / 媒体），因此优先展示。
func (s *Site) aggregateLocked() (State, Reason) {
	switch {
	case s.serverDegraded:
		return StateDegraded, ReasonProbeServer
	case s.houseDegraded:
		return StateDegraded, ReasonProbeHouse
	case s.write:
		return StateDegraded, ReasonWrite
	default:
		return StateNormal, ""
	}
}

// stop 停止探测任务并关闭全部订阅；可重复调用。
func (s *Site) stop() {
	s.stopOnce.Do(func() {
		s.scheduler.Stop()
		s.mu.Lock()
		s.closed = true
		for id, updates := range s.subscriptions {
			close(updates)
			delete(s.subscriptions, id)
		}
		s.mu.Unlock()
	})
}
