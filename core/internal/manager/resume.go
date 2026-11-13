package manager

import (
	"fmt"
	"os"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager/progress"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// loadLocalRebuild 检查本地是否已存在可续传的 run-*.swanlab 并一次扫描重建状态；
// 不存在或启用 SkipStore 时返回 (nil, nil)，已包含 FinishRecord 则拒绝续传。
func loadLocalRebuild(cfg config.Settings, run *runctx.RunContext) (*metrics.RebuildResult, error) {
	if cfg.SkipStore {
		return nil, nil
	}
	runFile, err := run.RunFile()
	if err != nil {
		return nil, fmt.Errorf("start: derive run file: %w", err)
	}
	info, err := os.Stat(runFile)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, nil
		}
		return nil, fmt.Errorf("start: stat run file %q: %w", runFile, err)
	}
	if info.IsDir() {
		return nil, fmt.Errorf("%w: run path %q is a directory", ErrInvalidStartRecord, runFile)
	}
	res, err := metrics.Rebuild(runFile)
	if err != nil {
		return nil, fmt.Errorf("%w: rebuild local store %q: %w", ErrInvalidStartRecord, runFile, err)
	}
	if res.Finished {
		return nil, fmt.Errorf("%w: run %q is already finished, use `swanlab sync` instead", ErrInvalidStartRecord, run.RunID)
	}
	return res, nil
}

// resolveResumeCheckpoint 读取并校验断点文件是否与当前可信本地副本范围一致；
// 缺失、损坏或越界时回退到文件头 (0, HeaderLen, false)。
func resolveResumeCheckpoint(
	checkpointFile, runID, runFile string,
	local *metrics.RebuildResult,
) (seq uint64, offset int64, valid bool) {
	fallback := func() (uint64, int64, bool) {
		return 0, store.HeaderLen, false
	}
	if local == nil || checkpointFile == "" {
		return fallback()
	}
	cp, ok, err := store.LoadCheckpoint(checkpointFile, runID)
	if err != nil {
		console.Debugf("resume: ignoring unusable checkpoint %s: %v", checkpointFile, err)
		return fallback()
	}
	if !ok {
		console.Debugf("resume: no matching checkpoint at %s, replaying from file head", checkpointFile)
		return fallback()
	}
	if cp.StoreSeq > local.Records || cp.Offset < store.HeaderLen || cp.Offset > local.EndOffset {
		console.Debugf("resume: checkpoint out of bounds (seq=%d/%d off=%d/%d), replaying from file head",
			cp.StoreSeq, local.Records, cp.Offset, local.EndOffset)
		return fallback()
	}
	if cp.StoreSeq == 0 && cp.Offset != store.HeaderLen {
		console.Debugf("resume: invalid zero-seq checkpoint offset %d, replaying from file head", cp.Offset)
		return fallback()
	}
	if cp.StoreSeq == local.Records && cp.Offset != local.EndOffset {
		console.Debugf("resume: checkpoint offset %d != end offset %d at seq %d, replaying from file head",
			cp.Offset, local.EndOffset, cp.StoreSeq)
		return fallback()
	}
	if cp.StoreSeq > 0 && cp.StoreSeq < local.Records {
		if cp.Offset >= local.EndOffset {
			console.Debugf("resume: checkpoint offset %d >= end offset %d with trailing records, replaying from file head",
				cp.Offset, local.EndOffset)
			return fallback()
		}
		if err = probeRecordBoundary(runFile, cp.Offset, local.EndOffset); err != nil {
			console.Debugf("resume: checkpoint offset %d is not a valid record boundary: %v", cp.Offset, err)
			return fallback()
		}
	}
	return cp.StoreSeq, cp.Offset, true
}

func probeRecordBoundary(runFile string, offset, limit int64) error {
	r, err := store.OpenOffset(runFile, offset, limit)
	if err != nil {
		return err
	}
	defer func() { _ = r.Close() }()

	payload, err := r.Scan()
	if err != nil {
		return err
	}
	rec := &recordv1.Record{}
	return proto.Unmarshal(payload, rec)
}

func countUploadable(runFile string, startOffset, endOffset int64) int64 {
	if startOffset >= endOffset {
		return 0
	}
	r, err := store.OpenOffset(runFile, startOffset, endOffset)
	if err != nil {
		return 0
	}
	defer func() { _ = r.Close() }()

	var n int64
	rec := &recordv1.Record{}
	for {
		payload, scanErr := r.Scan()
		if scanErr != nil {
			return n
		}
		rec.Reset()
		if proto.Unmarshal(payload, rec) == nil && record.IsImmediateUploadable(rec) {
			n++
		}
	}
}

func (s *Session) openStore(init sessionInit, resumeSeq uint64) (string, error) {
	if init.cfg.SkipStore {
		return "", nil
	}
	runFile, err := init.run.RunFile()
	if err != nil {
		return "", fmt.Errorf("derive run file: %w", err)
	}
	checkpointFile, err := init.run.CheckpointFile()
	if err != nil {
		return "", fmt.Errorf("derive checkpoint file: %w", err)
	}
	var writer *store.Writer
	if init.resume != nil {
		writer, err = store.OpenAppend(runFile, init.resume.EndOffset)
		if err != nil {
			return "", fmt.Errorf("open run file for append: %w", err)
		}
	} else {
		if err = os.MkdirAll(init.run.RunDir, 0o755); err != nil {
			return "", fmt.Errorf("create run dir: %w", err)
		}
		writer, err = store.Create(runFile)
		if err != nil {
			return "", fmt.Errorf("create run file: %w", err)
		}
		s.runFileCreated = true
	}
	s.runFile = runFile

	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{
		Queue:     queueSize(init.cfg.RecordBatch),
		ResumeSeq: resumeSeq,
		OnError: func(err error) {
			s.storeDegraded.Store(true)
			console.Warningf("store degraded, records will be uploaded without a local copy: %v", err)
			// 失去本地副本后不能再暂停，唤醒监听协程恢复投递。
			s.wakeHealth()
		},
	})
	if err != nil {
		_ = writer.Close()
		return "", fmt.Errorf("start store writer: %w", err)
	}
	s.store = async
	return checkpointFile, nil
}

func (s *Session) initOnlineProgress(init sessionInit, checkpointFile string) (
	onStored func(*recordv1.Record, uint64),
	replayFrom uint64,
	replayStartOffset int64,
	baseSeq uint64,
	cpValid bool,
) {
	replayStartOffset = store.HeaderLen
	if s.store == nil || init.cfg.NormalizedMode() != config.ModeOnline {
		return nil, 0, replayStartOffset, 0, false
	}
	if init.resume != nil {
		baseSeq, replayStartOffset, cpValid = resolveResumeCheckpoint(
			checkpointFile, init.run.RunID, s.runFile, init.resume,
		)
	}
	s.progress = progress.NewTracker(baseSeq)
	if init.resume != nil && baseSeq < init.resume.Records {
		replayFrom = baseSeq + 1
		s.backfillState = backfillState{
			pauseSeq: replayFrom,
			cursor:   backfillCursor{seq: baseSeq, offset: replayStartOffset},
		}
	}
	return s.progress.OnStored, replayFrom, replayStartOffset, baseSeq, cpValid
}
