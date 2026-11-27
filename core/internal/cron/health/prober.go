package health

import (
	"context"
	"sync"
	"sync/atomic"
	"time"
)

// ProbeFunc 返回服务是否可达；ctx 携带单次探测超时，实现必须尊重取消。
type ProbeFunc func(ctx context.Context) bool

// prober 按连续失败 / 成功阈值翻转一个服务的稳定状态。
// mu 串行化外部降级与 tick 的计数、状态翻转及回调；degraded 供并发读取。
type prober struct {
	probe        ProbeFunc
	timeout      time.Duration
	degradeAfter int
	recoverAfter int
	onFlip       func(degraded bool)

	mu       sync.Mutex
	degraded atomic.Bool

	failures  int
	successes int
}

type proberOptions struct {
	timeout      time.Duration
	degradeAfter int
	recoverAfter int
}

func newProber(probe ProbeFunc, opts proberOptions, onFlip func(degraded bool)) *prober {
	return &prober{
		probe:        probe,
		timeout:      opts.timeout,
		degradeAfter: opts.degradeAfter,
		recoverAfter: opts.recoverAfter,
		onFlip:       onFlip,
	}
}

// degradeNow 立即把服务标记为降级，供启动预检与共享站点合并使用。
// 同时清空失败与恢复进度，并与 tick 的状态结算互斥。
func (p *prober) degradeNow() {
	p.mu.Lock()
	defer p.mu.Unlock()
	p.failures = 0
	p.successes = 0
	p.degraded.Store(true)
	p.onFlip(true)
}

// tick 执行一次探测并按阈值更新稳定状态；由 cron 任务协程串行调用。
func (p *prober) tick(parent context.Context) {
	ctx, cancel := context.WithTimeout(parent, p.timeout)
	ok := p.probe(ctx)
	cancel()

	p.mu.Lock()
	defer p.mu.Unlock()
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
		// 降级期间继续探测，等待恢复阈值。
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
