// 本文件承载会话资源的装配、回滚与本地失效补传；
// 由 per-session 的 recovery 泵周期调用 intake.DrainRecovery。
package manager

import (
	"context"
	"errors"
	"fmt"
	"os"
	"time"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager/heartbeat"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	"github.com/swanhubx/swanlab/core/internal/uploader/sender"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// open 装配 store/intake/uploader/watcher/heartbeat；失败时由调用方执行 rollback。
func (s *Session) open(init sessionInit) error {
	var resumeSeq uint64
	var endSaves []*recordv1.Record
	if init.resume != nil {
		resumeSeq = init.resume.Records
		endSaves = init.resume.EndSaves
	}

	if err := s.openStore(init, resumeSeq); err != nil {
		return err
	}

	s.tracker = uploader.NewTracker()
	b := budget.New(budget.Limits{
		Records: init.cfg.UploadBudgetRecords,
		Bytes:   init.cfg.UploadBudgetBytes,
	})
	s.budget = b

	var storeIn chan<- store.Request
	var flushed func() uint64
	if s.store != nil {
		storeIn = s.store.In()
		flushed = func() uint64 {
			seq, _ := s.store.Flushed()
			return seq
		}
	}
	onStored, replayFrom, replayStartOffset := s.initOnlineProgress(init)
	var pendingDeferred int
	if len(endSaves) > 0 {
		pendingDeferred = buffer.DefaultPendingDeferred + len(endSaves)
	}
	intake, err := buffer.New(buffer.Options{
		SkipStore:       init.cfg.SkipStore,
		StoreOnly:       init.cfg.NormalizedMode() != config.ModeOnline,
		Store:           storeIn,
		StartSeq:        resumeSeq,
		PendingDeferred: pendingDeferred,
		UploadQueue:     queueSize(init.cfg.RecordBatch),
		Budget:          b,
		Degraded:        &s.storeDegraded,
		Flushed:         flushed,
		OnAdmit:         s.tracker.AddTotalRecords,
		OnPause:         func(uint64) { s.wakeBackfill() },
		OnStored:        onStored,
	})
	if err != nil {
		return fmt.Errorf("create intake: %w", err)
	}
	s.intake = intake
	s.intake.RestoreDeferred(endSaves)
	if replayFrom > 0 {
		s.intake.PauseAt(replayFrom)
	}

	if init.cfg.NormalizedMode() == config.ModeOnline {
		senderOpts := sender.Options{
			API:          init.api,
			Username:     init.workspace,
			Project:      init.project,
			ProjectID:    init.projectID,
			ExperimentID: init.experimentID,
			Files: sender.FileOptions{
				MediaDir:  init.run.MediaDir(),
				FilesDir:  init.run.FilesDir(),
				SkipStore: init.cfg.SkipStore,
				SaveSize:  init.cfg.SaveSize,
				SaveSplit: init.cfg.SaveSplit,
				SavePart:  init.cfg.SavePart,
				SaveBatch: int(init.cfg.SaveBatch),
				Tracker:   s.tracker,
				Putter:    init.http,
			},
		}
		mediaSender, err := sender.NewMedia(senderOpts)
		if err != nil {
			return fmt.Errorf("create media sender: %w", err)
		}
		saveSender, err := sender.NewSave(senderOpts)
		if err != nil {
			return fmt.Errorf("create save sender: %w", err)
		}

		up, err := uploader.New(uploader.Options{
			In:     intake.UploadOut(),
			Budget: b,
			Senders: map[record.Type]uploader.Sender{
				record.TypeColumn:    sender.NewColumn(senderOpts),
				record.TypeScalar:    sender.NewScalar(senderOpts),
				record.TypeLog:       sender.NewLog(senderOpts),
				record.TypeMedia:     mediaSender,
				record.TypeHistogram: sender.NewHistogram(senderOpts),
				record.TypeSave:      saveSender,
			},
			Batch:              int(init.cfg.RecordBatch),
			Interval:           time.Duration(init.cfg.RecordInterval * float64(time.Second)),
			Tracker:            s.tracker,
			Context:            s.ctx,
			OnRetryableFailure: s.onWriteFailure,
			OnDispatchSuccess:  s.onWriteSuccess,
			OnConfirm:          s.progress.Resolve,
			OnDrop:             s.progress.ResolveOne,
			OnDrained: func() {
				s.stopCheckpoint()
				s.flushCheckpoint()
				s.markDrained()
			},
		})
		if err != nil {
			return fmt.Errorf("create uploader: %w", err)
		}
		s.uploader = up
		s.drainFn = up.Wait
		s.releaseFn = func(ctx context.Context) error {
			return init.api.StopExperiment(ctx, api.StopExperimentParams{
				Username:     init.workspace,
				Project:      init.project,
				ExperimentID: init.experimentID,
				State:        s.stopState(),
			})
		}
	} else {
		// local/offline 没有上传阶段；收尾只需关闭本地 store 与 intake。
		s.drainFn = func(context.Context) error {
			s.markDrained()
			return nil
		}
	}

	if init.cfg.SkipStore {
		s.watcher = watcher.NullWatcher{}
	} else {
		fileWatcher, err := watcher.New(init.run.FilesDir(), s.onFileChanged)
		if err != nil {
			return fmt.Errorf("create watcher: %w", err)
		}
		s.watcher = fileWatcher
	}

	s.tracker.SetState(operationv1.CoreState_CORE_STATE_RUNNING)
	if init.cfg.NormalizedMode() == config.ModeOnline {
		// 心跳与站点探针只对云端有意义；local/offline/disabled 没有 api 客户端。
		s.heartbeat = heartbeat.New(init.experimentID, init.api.SendHeartbeat)
		s.heartbeat.Start(s.ctx)
		s.attachHealth(init.healthReg, init.healthOpts)
	}
	if s.store != nil {
		// skip_store 不会有本地写失败，无需补传协程。
		s.startRecovery()
	}
	if s.store != nil && s.uploader != nil {
		if init.resume != nil {
			if totalReplay := countUploadable(s.runFile, replayStartOffset, init.resume.EndOffset) + int64(len(endSaves)); totalReplay > 0 {
				s.tracker.AddTotalRecords(totalReplay)
			}
		}
		s.startBackfill()
		if replayFrom > 0 {
			s.wakeBackfill()
		}
	}
	s.startCheckpoint()
	return nil
}

// rollback 释放装配期已创建的资源；所有字段都允许为 nil。
func (s *Session) rollback() {
	if s.cancel != nil {
		s.cancel()
	}
	if s.heartbeat != nil {
		// s.cancel 已中断心跳的父 ctx；这里限制等待时间，避免慢心跳拖住回滚。
		ctx, cancel := context.WithTimeout(context.Background(), rollbackTimeout)
		s.heartbeat.Stop(ctx)
		cancel()
	}
	if s.watcher != nil {
		s.watcher.Stop()
	}
	s.stopRecovery()
	s.stopBackfill()
	s.stopCheckpoint()
	if s.uploader != nil {
		_ = s.uploader.Close()
	}
	if s.intake != nil {
		ctx, cancel := context.WithTimeout(context.Background(), rollbackTimeout)
		_ = s.intake.Close(ctx)
		cancel()
	}
	if s.store != nil {
		_ = s.store.Close()
	}
	if s.runFileCreated {
		_ = os.Remove(s.runFile)
		if s.checkpointFile != "" {
			_ = os.Remove(s.checkpointFile)
		}
	}
	s.closeHealth()
}

// queueSize 是上传通道与 store 通道容量：单批记录最多因自动列翻倍。
func queueSize(batch int32) int { return queueFactor * int(batch) }

func firstNonEmpty(values ...string) string {
	for _, value := range values {
		if value != "" {
			return value
		}
	}
	return ""
}

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
