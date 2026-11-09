package health

import (
	"context"
	"sync/atomic"
	"time"
)

// ProbeFunc 返回平面是否可达；ctx 携带单次探测超时，实现必须尊重取消。
type ProbeFunc func(ctx context.Context) bool

// prober 按连续失败 / 成功阈值翻转一个平面的稳定状态。
//
// 计数器只在探测协程内读写；reachable / degraded 供 Snapshot 并发读取。
type prober struct {
	probe        ProbeFunc
	timeout      time.Duration
	degradeAfter int
	recoverAfter int
	onFlip       func(degraded bool)

	reachable atomic.Bool
	degraded  atomic.Bool

	failures  int
	successes int
}

type proberOptions struct {
	timeout      time.Duration
	degradeAfter int
	recoverAfter int
}

func newProber(probe ProbeFunc, opts proberOptions, onFlip func(degraded bool)) *prober {
	p := &prober{
		probe:        probe,
		timeout:      opts.timeout,
		degradeAfter: opts.degradeAfter,
		recoverAfter: opts.recoverAfter,
		onFlip:       onFlip,
	}
	// 首轮探测前保持乐观，避免新站点在第一个周期内被误判为降级。
	p.reachable.Store(true)
	return p
}

// seedDegraded 写入启动预检结论；只允许在探测协程启动前调用。
func (p *prober) seedDegraded() {
	p.reachable.Store(false)
	p.degraded.Store(true)
}

// tick 执行一次探测并按阈值更新稳定状态；由 cron 任务协程串行调用。
func (p *prober) tick(parent context.Context) {
	ctx, cancel := context.WithTimeout(parent, p.timeout)
	ok := p.probe(ctx)
	cancel()

	p.reachable.Store(ok)
	switch {
	case ok && p.degraded.Load():
		p.failures = 0
		p.successes++
		if p.successes >= p.recoverAfter {
			p.successes = 0
			p.degraded.Store(false)
			p.onFlip(false)
		}
	case ok:
		p.failures = 0
	case p.degraded.Load():
		// 降级期间继续探测，成功计数已在上面处理。
		p.successes = 0
	default:
		p.successes = 0
		p.failures++
		if p.failures >= p.degradeAfter {
			p.failures = 0
			p.degraded.Store(true)
			p.onFlip(true)
		}
	}
}
