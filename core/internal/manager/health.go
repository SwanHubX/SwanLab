// 本文件承载会话与站点健康状态的接线：探针注册、暂停/恢复与降级信号。
package manager

import (
	"context"

	"github.com/swanhubx/swanlab/core/internal/cron/health"
)

// writeDegradeRounds 是连续可重试失败达到降级所需的轮数。
const writeDegradeRounds = 3

// attachHealth 获取站点引用、独立订阅并启动健康监听；启动预检已降级的站点在返回前完成暂停。
func (s *Session) attachHealth(reg *health.Registry, opts health.SiteOptions) {
	if reg == nil {
		return
	}
	s.healthReg = reg
	s.site = reg.Acquire(opts)
	s.healthUpdates, s.healthUnsubscribe = s.site.Subscribe()

	ctx, cancel := context.WithCancel(context.Background())
	s.healthStop = cancel
	// 订阅先于初始结算，启动监听前发生的变化也会保留在独立 channel 中。
	s.reconcileHealth()

	s.healthWG.Add(1)
	go func() {
		defer s.healthWG.Done()
		s.watchHealth(ctx)
	}()
}

// watchHealth 串行处理站点状态变化与 store 降级事件，避免并发决定暂停/恢复。
func (s *Session) watchHealth(ctx context.Context) {
	for {
		select {
		case <-ctx.Done():
			return
		case <-s.healthWake:
		case _, ok := <-s.healthUpdates:
			if !ok {
				return
			}
		}
		s.reconcileHealth()
	}
}

// reconcileHealth 把站点状态映射为上传暂停/恢复，并驱动补投启停。
// 写入侧降级（reason=write）保持投递，否则暂停后永远等不到成功信号。
func (s *Session) reconcileHealth() {
	if s.site == nil || s.uploader == nil {
		return
	}
	if s.shouldSuspendUploads() {
		s.uploader.Suspend()
		s.backfill.CancelAttempt()
		return
	}
	s.uploader.Resume()
	// 健康状态变化后让补投按最新状态重新开始。
	s.backfill.CancelAttempt()
	s.backfill.Wake()
}

// shouldSuspendUploads 判断当前是否必须暂停上传投递。
func (s *Session) shouldSuspendUploads() bool {
	if s.site == nil || s.cfg.SkipStore || s.storeDegraded.Load() {
		return false
	}
	snap := s.site.Snapshot()
	return snap.State == health.StateDegraded && snap.Reason != health.ReasonWrite
}

// wakeHealth 唤醒监听协程重新评估暂停状态；监听启动前的事件由初始 reconcile 覆盖。
func (s *Session) wakeHealth() {
	select {
	case s.healthWake <- struct{}{}:
	default:
	}
}

// onWriteFailure 在连续可重试失败达到阈值时注入写入侧降级。
func (s *Session) onWriteFailure(rounds int) {
	if rounds >= writeDegradeRounds && s.site != nil {
		s.site.Report(health.ReasonWrite, true)
	}
}

// onWriteSuccess 在成功调度一轮后解除写入侧降级。
func (s *Session) onWriteSuccess() {
	if s.site != nil {
		s.site.Report(health.ReasonWrite, false)
	}
}

// closeHealth 幂等停止监听、注销并释放站点引用；必须在 uploader 关闭后调用。
func (s *Session) closeHealth() {
	s.healthOnce.Do(func() {
		if s.healthStop != nil {
			s.healthStop()
		}
		s.healthWG.Wait()
		if s.healthUnsubscribe != nil {
			s.healthUnsubscribe()
		}
		if s.healthReg != nil && s.site != nil {
			s.healthReg.Release(s.site)
		}
	})
}
