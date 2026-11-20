// Package recovery 承载本地存储降级后的回收队列补传泵，周期调用 DrainRecovery，
// 并在会话收尾时尽最大努力排空回收队列。
package recovery

import (
	"context"
	"errors"
	"sync"
	"sync/atomic"
	"time"

	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

const (
	// defaultInterval 是本地失效后重试补传的周期；降级未发生时只空转。
	defaultInterval = 500 * time.Millisecond
	// defaultPollInterval 是收尾阶段等待回收队列腾出上传槽位的轮询间隔。
	defaultPollInterval = 10 * time.Millisecond
)

// Source 抽象回收队列的排空能力（生产环境由 *buffer.Intake 实现）。
type Source interface {
	DrainRecovery(ctx context.Context) (remaining int, err error)
}

// Options 配置本地失效补传泵。
type Options struct {
	Ctx          context.Context
	Source       Source
	Degraded     *atomic.Bool
	Interval     time.Duration
	PollInterval time.Duration
}

// Pump 管理 per-session 的本地失效补传协程与停机排空；nil 接收者上的方法均为安全空操作。
type Pump struct {
	ctx          context.Context
	source       Source
	degraded     *atomic.Bool
	interval     time.Duration
	pollInterval time.Duration

	stop      chan struct{}
	done      chan struct{}
	startOnce sync.Once
	stopOnce  sync.Once
}

// New 创建补传泵；需调用 Start 启动后台协程。
func New(opts Options) *Pump {
	ctx := opts.Ctx
	if ctx == nil {
		ctx = context.Background()
	}
	interval := opts.Interval
	if interval <= 0 {
		interval = defaultInterval
	}
	pollInterval := opts.PollInterval
	if pollInterval <= 0 {
		pollInterval = defaultPollInterval
	}
	return &Pump{
		ctx:          ctx,
		source:       opts.Source,
		degraded:     opts.Degraded,
		interval:     interval,
		pollInterval: pollInterval,
	}
}

// Start 启动后台补传协程；只能启动一次，且在 Intake.Close 之前必须调用 Stop。
func (p *Pump) Start() {
	if p == nil || p.source == nil {
		return
	}
	p.startOnce.Do(func() {
		p.stop = make(chan struct{})
		p.done = make(chan struct{})
		go p.loop()
	})
}

// Stop 停止并等待补传协程退出；未启动或已停止时可重复安全调用。
func (p *Pump) Stop() {
	if p == nil {
		return
	}
	p.stopOnce.Do(func() {
		if p.stop != nil {
			close(p.stop)
		}
	})
	if p.done != nil {
		<-p.done
	}
}

// Drain 在收尾阶段尽量排空回收队列；队列仍满时按 poll 间隔重试，
// 直到队列为空、返回错误或 ctx 到期。
func (p *Pump) Drain(ctx context.Context) {
	if p == nil || p.source == nil {
		return
	}
	for {
		remaining, err := p.source.DrainRecovery(ctx)
		if remaining == 0 || err != nil {
			return
		}
		timer := time.NewTimer(p.pollInterval)
		select {
		case <-ctx.Done():
			timer.Stop()
			return
		case <-timer.C:
		}
	}
}

func (p *Pump) loop() {
	defer close(p.done)
	ticker := time.NewTicker(p.interval)
	defer ticker.Stop()
	for {
		select {
		case <-p.ctx.Done():
			return
		case <-p.stop:
			return
		case <-ticker.C:
		}
		if p.degraded != nil && !p.degraded.Load() {
			continue
		}
		remaining, err := p.source.DrainRecovery(p.ctx)
		switch {
		case errors.Is(err, buffer.ErrClosed), errors.Is(err, context.Canceled):
			return
		case err != nil:
			console.Warningf("recovery: retry local store fallback: %v", err)
		case remaining > 0:
			console.Debugf("recovery: %d record(s) waiting for upload capacity", remaining)
		}
	}
}
