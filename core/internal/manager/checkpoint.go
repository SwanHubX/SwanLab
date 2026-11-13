package manager

import (
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/store"
)

const defaultCheckpointInterval = time.Second

// startCheckpoint 启动周期持久化上传断点的后台协程；skip_store 或无本地库时跳过。
func (s *Session) startCheckpoint() {
	if s.checkpointFile == "" || s.store == nil || s.run == nil {
		return
	}
	interval := s.checkpointInterval
	if interval <= 0 {
		interval = defaultCheckpointInterval
	}
	var ctxDone <-chan struct{}
	if s.ctx != nil {
		ctxDone = s.ctx.Done()
	}
	s.checkpointStop = make(chan struct{})
	s.checkpointDone = make(chan struct{})
	go func() {
		defer close(s.checkpointDone)
		ticker := time.NewTicker(interval)
		defer ticker.Stop()
		for {
			select {
			case <-ctxDone:
				return
			case <-s.checkpointStop:
				return
			case <-ticker.C:
				s.persistCheckpoint()
			}
		}
	}()
}

// stopCheckpoint 停止断点协程并等待其退出；未启动时为安全空操作，可重复调用。
func (s *Session) stopCheckpoint() {
	if s.checkpointStop == nil {
		return
	}
	s.checkpointOnce.Do(func() {
		close(s.checkpointStop)
	})
	<-s.checkpointDone
}

// flushCheckpoint 在停机或排空完成后刷写一次最终断点。
func (s *Session) flushCheckpoint() {
	s.persistCheckpoint()
}

// persistCheckpoint 把当前已连续确认且已落盘的进度写入 checkpoint.json。
// 在线模式取 min(progress.watermark(), flushedSeq)；本地/离线模式直接取 flushedSeq。
// 写入失败仅告警一次并留待后续周期重试，不中断会话运行。
func (s *Session) persistCheckpoint() {
	if s.checkpointFile == "" || s.store == nil || s.run == nil {
		return
	}
	s.checkpointMu.Lock()
	defer s.checkpointMu.Unlock()

	flushedSeq, flushedOff := s.store.Flushed()
	targetSeq := flushedSeq
	if s.progress != nil {
		if w := s.progress.Watermark(); w < targetSeq {
			targetSeq = w
		}
	}
	if targetSeq == 0 {
		return
	}
	if s.checkpointHasWrite && targetSeq == s.checkpointLastSeq {
		return
	}

	var (
		offset int64
		err    error
	)
	if targetSeq == flushedSeq {
		offset = flushedOff
	} else {
		offset, err = store.EndOffset(s.runFile, s.store.Index(), targetSeq, flushedOff)
		if err != nil {
			s.warnCheckpoint(err)
			return
		}
	}

	cp := store.Checkpoint{
		Version:   store.CheckpointVersion,
		RunID:     s.run.RunID,
		StoreSeq:  targetSeq,
		Offset:    offset,
		UpdatedAt: time.Now(),
	}
	if err = store.WriteCheckpoint(s.checkpointFile, cp); err != nil {
		s.warnCheckpoint(err)
		return
	}
	s.checkpointWarned = false
	s.checkpointLastSeq = targetSeq
	s.checkpointHasWrite = true
	if s.checkpointTicked != nil {
		s.checkpointTicked()
	}
}

func (s *Session) warnCheckpoint(err error) {
	if s.checkpointWarned {
		return
	}
	s.checkpointWarned = true
	console.Warningf("checkpoint: write %s: %v", s.checkpointFile, err)
}
