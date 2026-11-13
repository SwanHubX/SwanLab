// 本文件承载会话资源的装配、回滚与本地失效补传；
// 由 per-session 的 recovery 泵周期调用 intake.DrainRecovery。
package manager

import (
	"context"
	"fmt"
	"os"
	"time"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager/backfill"
	"github.com/swanhubx/swanlab/core/internal/manager/checkpoint"
	"github.com/swanhubx/swanlab/core/internal/manager/heartbeat"
	"github.com/swanhubx/swanlab/core/internal/manager/recovery"
	"github.com/swanhubx/swanlab/core/internal/manager/resume"
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

	opened, err := resume.OpenStore(resume.StoreOptions{
		SkipStore: init.cfg.SkipStore,
		Run:       init.run,
		Resume:    init.resume,
		Queue:     queueSize(init.cfg.RecordBatch),
		OnError: func(err error) {
			s.storeDegraded.Store(true)
			console.Warningf("store degraded, records will be uploaded without a local copy: %v", err)
			// 失去本地副本后不能再暂停，唤醒监听协程恢复投递。
			s.wakeHealth()
		},
	})
	if err != nil {
		return err
	}
	s.store = opened.Writer
	s.runFile = opened.RunFile
	s.runFileCreated = opened.Created

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
	prog := resume.InitOnlineProgress(resume.ProgressOptions{
		Mode:           init.cfg.NormalizedMode(),
		RunID:          init.run.RunID,
		RunFile:        s.runFile,
		CheckpointFile: opened.CheckpointFile,
		HasStore:       s.store != nil,
		Resume:         init.resume,
	})
	s.progress = prog.Tracker
	if init.cfg.NormalizedMode() == config.ModeOnline && s.progress != nil {
		s.checkpoint = checkpoint.New(checkpoint.Options{
			Ctx:            s.ctx,
			File:           opened.CheckpointFile,
			RunID:          init.run.RunID,
			RunFile:        s.runFile,
			Store:          s.store,
			Watermark:      s.progress.Watermark,
			InitialSeq:     prog.CheckpointSeq,
			InitialWritten: prog.CheckpointValid && prog.CheckpointSeq > 0,
		})
	}
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
		OnPause:         func(uint64) { s.backfill.Wake() },
		OnStored:        prog.OnStored,
	})
	if err != nil {
		return fmt.Errorf("create intake: %w", err)
	}
	s.intake = intake
	s.intake.RestoreDeferred(endSaves)
	if prog.ReplayFrom > 0 {
		s.intake.PauseAt(prog.ReplayFrom)
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
				s.checkpoint.Stop()
				s.checkpoint.Flush()
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
		s.recovery = recovery.New(recovery.Options{
			Ctx:      s.ctx,
			Source:   s.intake,
			Degraded: &s.storeDegraded,
		})
		s.recovery.Start()
	}
	if s.store != nil && s.uploader != nil {
		if init.resume != nil {
			if totalReplay := resume.CountUploadable(s.runFile, prog.ReplayStartOffset, init.resume.EndOffset) + int64(len(endSaves)); totalReplay > 0 {
				s.tracker.AddTotalRecords(totalReplay)
			}
		}
		var initCursor backfill.Cursor
		if prog.ReplayFrom > 0 {
			initCursor = backfill.Cursor{Seq: prog.ReplayFrom - 1, Offset: prog.ReplayStartOffset}
		}
		var onBackfillStored func(*recordv1.Record, uint64)
		if s.progress != nil {
			onBackfillStored = s.progress.OnStored
		}
		s.backfill = backfill.New(backfill.Options{
			Ctx:           s.ctx,
			Intake:        s.intake,
			Store:         s.store,
			RunFile:       s.runFile,
			BatchSize:     s.cfg.RecordBatch,
			SkipStore:     s.cfg.SkipStore,
			Degraded:      &s.storeDegraded,
			IsRunning:     func() bool { return s.stateSnapshot() == RunStateRunning },
			IsSuspended:   s.shouldSuspendUploads,
			OnStored:      onBackfillStored,
			InitialPause:  prog.ReplayFrom,
			InitialCursor: initCursor,
		})
		s.backfill.Start()
		if prog.ReplayFrom > 0 {
			s.backfill.Wake()
		}
	}
	s.checkpoint.Start()
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
	s.recovery.Stop()
	s.backfill.Stop()
	s.checkpoint.Stop()
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
		s.checkpoint.RemoveFile()
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
