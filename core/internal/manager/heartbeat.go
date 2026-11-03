package manager

import (
	"context"
	"sync"
	"sync/atomic"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// heartbeatInterval 是会话心跳周期，与 Python 侧保持一致。
const heartbeatInterval = 10 * time.Minute

type heartbeatSender func(context.Context, string) error

// heartbeat 维持会话活跃状态：启动时立即发送一次，之后按固定周期重发。
// 发送失败只告警，不影响数据主链路。
type heartbeat struct {
	experimentID string
	send         heartbeatSender
	interval     time.Duration

	started   atomic.Bool
	startOnce sync.Once
	stopOnce  sync.Once
	stopCh    chan struct{}
	doneCh    chan struct{}
}

func newHeartbeat(experimentID string, send heartbeatSender) *heartbeat {
	return &heartbeat{
		experimentID: experimentID,
		send:         send,
		interval:     heartbeatInterval,
		stopCh:       make(chan struct{}),
		doneCh:       make(chan struct{}),
	}
}

// start 启动心跳协程；重复调用无副作用。
func (h *heartbeat) start(ctx context.Context) {
	h.startOnce.Do(func() {
		h.started.Store(true)
		go h.loop(ctx)
	})
}

// stop 停止心跳并等待协程退出；可重复调用，未启动时不会阻塞。
func (h *heartbeat) stop() {
	h.stopOnce.Do(func() {
		close(h.stopCh)
	})
	if h.started.Load() {
		<-h.doneCh
	}
}

func (h *heartbeat) loop(ctx context.Context) {
	defer close(h.doneCh)

	ticker := time.NewTicker(h.interval)
	defer ticker.Stop()

	if h.stopped(ctx) {
		return
	}
	h.beat(ctx)
	for {
		select {
		case <-h.stopCh:
			return
		case <-ctx.Done():
			return
		case <-ticker.C:
			h.beat(ctx)
		}
	}
}

func (h *heartbeat) stopped(ctx context.Context) bool {
	select {
	case <-h.stopCh:
		return true
	case <-ctx.Done():
		return true
	default:
		return false
	}
}

func (h *heartbeat) beat(ctx context.Context) {
	if ctx.Err() != nil {
		return
	}
	if err := h.send(ctx, h.experimentID); err != nil {
		console.Warningf("heartbeat failed: %v", err)
	}
}
