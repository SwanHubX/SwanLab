// 本文件承载 DeliverRunFinish 与进程退出的收尾序列。三条不变量：
//   - FinishRecord 必须先提交到本地副本，之后才允许关闭 store；
//   - finishCommitted 保证重试时不会重复写入 FinishRecord；
//   - pendingLog 非空时不得关闭 intake，否则崩溃 error_log 无法补传。
package manager

import (
	"context"
	"errors"
	"fmt"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

// finish 执行或续跑 DeliverRunFinish 的收尾序列：
//
//	DRAINING → watcher.Stop → 保存终态并落盘 FinishRecord → 阻塞刷盘 →
//	补传崩溃 error_log → Intake.Close（冲刷 END 暂存与 recovery）→ RequestFinish
//
// 收尾序列失败后会话进入 FAILED，但本地 store 与 intake 保持可用：重复投递会
// 从失败点续跑，Confirm 也会用新的 ctx 重试同样的步骤。FinishRecord 落盘与
// error_log 入队是释放前必须完成的步骤，不能因为暂时性失败被跳过。
func (s *Session) finish(ctx context.Context, rec *runv1.FinishRecord) error {
	s.finishMu.Lock()
	defer s.finishMu.Unlock()

	switch s.stateSnapshot() {
	case RunStateRunning:
		started, err := s.beginDrainOwned()
		if err != nil {
			return err
		}
		if !started {
			return nil
		}
		s.stopBackfill()
		// watcher.Stop 会等待在途回调并同步检查一次文件状态；此时状态已是 DRAINING，
		// 防抖窗口内的残留回调会被状态门拒绝并丢弃，与 Python 的 _store_finish 一致。
		if s.watcher != nil {
			s.watcher.Stop()
		}
		// 终态必须在 Reserve 之前保存：落盘失败时 stop_experiment 仍按真实终态上报。
		s.prepareFinish(rec)
	case RunStateFailed:
		// 上次收尾中断，续跑剩余步骤；终态记录已保存，忽略新的 FinishRecord。
	default:
		// DRAINING/DRAINED 也可能由 Shutdown 切换而来，此时没有终态记录。
		// 只有已提交过 FinishRecord，重复投递才按幂等成功返回。
		if !s.finishCommitted {
			return fmt.Errorf("%w: run is shutting down before finish was recorded", ErrRunFailed)
		}
		return nil
	}

	if err := s.completeFinish(ctx); err != nil {
		s.fail(err)
		return s.failure()
	}
	return nil
}

// finishCommittedNow 在 finishMu 下读取终态提交标记，与 DeliverRunFinish 和
// Shutdown 的收尾序列互斥，避免读到中间状态。
func (s *Session) finishCommittedNow() bool {
	s.finishMu.Lock()
	defer s.finishMu.Unlock()
	return s.finishCommitted
}

// prepareFinish 保存终态记录并构建崩溃 error_log。
//
// 崩溃 error_log 按 Python 的设计只上传不落盘；暂时入队失败时由 pendingLog
// 保留，交给后续 DeliverRunFinish/Confirm 重试，不能静默丢弃。
func (s *Session) prepareFinish(rec *runv1.FinishRecord) {
	finish := proto.Clone(rec).(*runv1.FinishRecord)

	s.mu.Lock()
	defer s.mu.Unlock()
	s.finishRecord = finish
	// local/offline 不做上传，崩溃日志只保留在 FinishRecord.error 中。
	if s.cfg.NormalizedMode() != config.ModeOnline {
		return
	}
	if finish.GetState() == runv1.RunState_RUN_STATE_FINISHED {
		return
	}
	line := finish.GetError()
	if line == "" {
		line = finishUnknownError
	}
	s.pendingLog = builder.BuildLogRecord(s.counter, s.epoch, &terminalv1.LogRecord{
		Line:      line,
		Level:     terminalv1.LogLevel_LOG_LEVEL_ERROR,
		Timestamp: finish.GetFinishedAt(),
	})
}

// completeFinish 依次完成终态记录提交、本地刷盘、暂停区间补投、崩溃日志入队与
// 通道关闭；各步骤均可重入，便于失败后重试。
func (s *Session) completeFinish(ctx context.Context) error {
	if !s.finishCommitted {
		if err := s.commitFinish(ctx); err != nil {
			return err
		}
	}

	// 本地刷盘不设超时，必须等 AsyncWriter 全部写入，且关闭要排在 FinishRecord 提交之后。
	// 刷盘失败时 online 由上传兜底，local/offline 返回失败。
	if s.store != nil {
		if err := s.store.Close(); err != nil {
			s.storeDegraded.Store(true)
			console.Warningf("finish: local store did not flush completely: %v", err)
		}
		if s.uploader == nil {
			s.stopCheckpoint()
			s.flushCheckpoint()
		}
	}

	// 本地记录已全部刷盘，补完暂停区间后才能关闭 intake。
	if err := s.completeBackfill(ctx); err != nil {
		return err
	}

	if err := s.commitCrashLog(ctx); err != nil {
		return err
	}
	s.stopRecovery()
	closeErr := s.closeIntake(ctx)
	if s.cfg.NormalizedMode() != config.ModeOnline && s.storeDegraded.Load() {
		// local/offline 无云端兜底：关闭 intake 后返回失败，不 markDrained。
		if closeErr != nil {
			console.Warningf("finish: close intake after local store failure: %v", closeErr)
		}
		return fmt.Errorf("%w: local copy is incomplete", ErrRunFailed)
	}
	if closeErr != nil {
		return closeErr
	}
	if s.uploader == nil {
		s.markDrained()
	}
	return nil
}

// commitFinish 把 FinishRecord 写入本地副本；失败时保留 finishCommitted=false，
// 由重复的 DeliverRunFinish 或 Confirm 用新的 ctx 重试。
func (s *Session) commitFinish(ctx context.Context) error {
	s.mu.Lock()
	finish := s.finishRecord
	s.mu.Unlock()
	if finish == nil {
		return fmt.Errorf("%w: finish record is missing", errFinishRecord)
	}
	if s.intake == nil {
		return fmt.Errorf("%w: session has no intake", errFinishRecord)
	}
	if s.storeDegraded.Load() {
		// online 仍可依赖上传兜底；local/offline 的本地副本是唯一事实源，不能假装成功。
		if s.cfg.NormalizedMode() != config.ModeOnline {
			return fmt.Errorf("%w: local store is degraded, finish record is not persisted", ErrRunFailed)
		}
		console.Warning("finish: local store is degraded, finish record is not kept locally")
		s.finishCommitted = true
		return nil
	}

	ticket, err := s.intake.Reserve(ctx, buffer.Demand{Store: 1})
	if err != nil {
		return fmt.Errorf("%w: reserve: %w", errFinishRecord, err)
	}
	if err = s.intake.Commit(ticket, []*recordv1.Record{builder.BuildFinishRecord(finish)}); err != nil {
		return fmt.Errorf("%w: commit: %w", errFinishRecord, err)
	}
	s.finishCommitted = true
	return nil
}

// commitCrashLog 把崩溃 error_log 投入上传通道。暂时失败（通道满、ctx 到期）时
// 保留记录供 Confirm 重试；只有确定性错误才告警丢弃，避免会话永远无法释放。
func (s *Session) commitCrashLog(ctx context.Context) error {
	if s.cfg.NormalizedMode() != config.ModeOnline {
		return nil
	}
	if s.pendingLog == nil || s.intake == nil {
		return nil
	}
	if s.storeDegraded.Load() {
		// 回收队列中是更早的记录，error_log 不能插到它们前面。
		s.drainRecovery(ctx)
	}

	err := s.intake.CommitUploadOnly(ctx, s.pendingLog)
	switch {
	case err == nil:
		s.pendingLog = nil
		return nil
	case errors.Is(err, budget.ErrTooLarge), errors.Is(err, buffer.ErrClosed):
		// 重试不会成功的确定性错误：告警后终结，不能让收尾永久卡住。
		console.Warningf("finish: dropping crash error log: %v", err)
		s.pendingLog = nil
		return nil
	default:
		return fmt.Errorf("upload crash error log: %w", err)
	}
}

// closeIntake 冲刷 END 暂存与 recovery、关闭上传通道，并通知 uploader 收尾；
// 只有 Close 成功才允许 RequestFinish，避免 uploader 在通道未封闭时退出。
func (s *Session) closeIntake(ctx context.Context) error {
	if s.intake == nil {
		return nil
	}
	if err := s.intake.Close(ctx); err != nil {
		return fmt.Errorf("%w: close intake: %w", ErrRunFailed, err)
	}
	if s.uploader != nil {
		s.uploader.RequestFinish()
	}
	return nil
}

// resumeFinish 在 Confirm 命中 FAILED 时续跑被中断的收尾序列：重试 FinishRecord
// 落盘、崩溃 error_log 入队与 Intake.Close，并等待 uploader 真正排空。
// 全部成功后才把状态恢复为 DRAINED；任一步仍失败时保留 FAILED 供下次重试，
// 调用方不得在这种情况下上报实验终态。
func (s *Session) resumeFinish(ctx context.Context) error {
	s.finishMu.Lock()
	defer s.finishMu.Unlock()

	if s.stateSnapshot() != RunStateFailed {
		return nil
	}
	if s.intake == nil {
		return s.failure()
	}
	if err := s.completeFinish(ctx); err != nil {
		s.fail(err)
		return s.failure()
	}
	if s.uploader != nil {
		if err := s.uploader.Wait(ctx); err != nil {
			return fmt.Errorf("confirm: wait for drain: %w", err)
		}
	}
	if !s.markRecovered() {
		return s.failure()
	}
	return nil
}

// stopState 返回上报给云端的实验终态；未记录 FinishRecord 时按正常结束处理。
func (s *Session) stopState() api.RunState {
	s.mu.Lock()
	defer s.mu.Unlock()
	return runStateFromProto(s.finishRecord.GetState())
}

// confirmMessage 是按本地副本完整性生成的 stop_experiment 失败文案。
func (s *Session) confirmMessage() string {
	switch {
	case s.cfg.SkipStore:
		return "Failed to finish run, and no local copy was kept (skip_store); the run state on the cloud is unconfirmed."
	case s.storeDegraded.Load():
		return "Failed to finish run, and the local copy is incomplete; the run state on the cloud is unconfirmed."
	default:
		return "Failed to finish run, but it has been saved locally."
	}
}

// runStateFromProto 把 proto 枚举转换为云端上报的终态。
// 与 Python 的 stop_experiment 一致：只有 crashed/aborted 是特殊值，
// 其余（含 RUNNING 与未知值）一律按正常结束处理。
func runStateFromProto(state runv1.RunState) api.RunState {
	switch state {
	case runv1.RunState_RUN_STATE_CRASHED:
		return api.RunStateCrashed
	case runv1.RunState_RUN_STATE_ABORTED:
		return api.RunStateAborted
	default:
		return api.RunStateFinished
	}
}

// shutdown 执行进程退出时的有界收尾：不写 FinishRecord，也不调用 stop_experiment。
// 顺序与 DeliverRunFinish 对齐，但全部等待受 ctx 约束，供 Manager.Shutdown 汇总。
// 持有 finishMu 与在途的 DeliverRunFinish/Confirm 互斥：退出钩子在 GracefulStop
// 之前运行，否则可能一边关闭本地通道一边写入。
func (s *Session) shutdown(ctx context.Context) {
	s.finishMu.Lock()
	defer s.finishMu.Unlock()

	_ = s.beginDrain()
	s.stopBackfill()
	if s.heartbeat != nil {
		s.heartbeat.stop(ctx)
	}
	if s.watcher != nil {
		s.watcher.Stop()
	}
	s.stopRecovery()
	if s.store != nil {
		if err := s.store.CloseContext(ctx); err != nil {
			console.Warningf("shutdown: close local store: %v", err)
		}
	}
	// 上次收尾中断遗留的崩溃 error_log 在封闭通道前补投，失败只能放弃。
	if err := s.commitCrashLog(ctx); err != nil {
		console.Warningf("shutdown: upload crash error log: %v", err)
	}
	closed := true
	if s.intake != nil {
		if err := s.intake.Close(ctx); err != nil {
			console.Warningf("shutdown: close intake: %v", err)
			closed = false
		}
	}
	if closed && s.uploader != nil {
		s.uploader.RequestFinish()
	}
	if s.uploader != nil {
		if err := s.uploader.Wait(ctx); err != nil {
			console.Warningf("shutdown: wait for drain: %v", err)
		}
	}
	s.stopCheckpoint()
	s.flushCheckpoint()
	if s.cancel != nil {
		s.cancel()
	}
	s.closeHealth()
}
