// 本文件承载本地存储失效后的补传协程：AsyncWriter 只置位 storeDegraded，
// 由 per-session 的 recovery 泵周期调用 intake.DrainRecovery。
package manager

import (
	"context"
	"errors"
	"time"

	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

const (
	// recoveryInterval 是本地失效后重试补传的周期；降级未发生时只空转。
	recoveryInterval = 500 * time.Millisecond
	// recoveryPollInterval 是收尾阶段等待回收队列腾出上传槽位的轮询间隔。
	recoveryPollInterval = 10 * time.Millisecond
)

// startRecovery 启动 per-session 的补传协程；只在 open 装配成功后调用一次。
func (s *Session) startRecovery() {
	s.recoveryStop = make(chan struct{})
	s.recoveryDone = make(chan struct{})
	go s.recoveryLoop()
}

// recoveryLoop 在 storeDegraded 为真时重试回收队列；Intake.Close 之前必须停止。
func (s *Session) recoveryLoop() {
	defer close(s.recoveryDone)
	ticker := time.NewTicker(recoveryInterval)
	defer ticker.Stop()
	for {
		select {
		case <-s.ctx.Done():
			return
		case <-s.recoveryStop:
			return
		case <-ticker.C:
		}
		if !s.storeDegraded.Load() {
			continue
		}
		remaining, err := s.intake.DrainRecovery(s.ctx)
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

// stopRecovery 停止并等待补传协程退出；可重复调用。
func (s *Session) stopRecovery() {
	s.recoveryOnce.Do(func() {
		if s.recoveryStop != nil {
			close(s.recoveryStop)
		}
	})
	if s.recoveryDone != nil {
		<-s.recoveryDone
	}
}

// drainRecovery 在收尾阶段尽量排空回收队列；队列仍满时按 poll 间隔重试，
// 直到队列为空或 ctx 到期（由调用方的 RPC deadline 兜底，避免永久阻塞）。
func (s *Session) drainRecovery(ctx context.Context) {
	for {
		remaining, err := s.intake.DrainRecovery(ctx)
		if remaining == 0 || err != nil {
			return
		}
		timer := time.NewTimer(recoveryPollInterval)
		select {
		case <-ctx.Done():
			timer.Stop()
			return
		case <-timer.C:
		}
	}
}
