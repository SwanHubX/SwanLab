// Package resume 负责崩溃续传前的本地 store 重建、checkpoint 校验、
// 追加写打开以及在线重放水位初始化。
package resume

import (
	"errors"
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

// ErrInvalidRunFile 表示已有的本地 run 文件不合法（如路径是目录、文件损坏或已完成运行）。
var ErrInvalidRunFile = errors.New("invalid local run file")

// LoadLocalRebuild 检查本地是否已存在可续传的 run-*.swanlab 并一次扫描重建状态；
// 不存在或启用 SkipStore 时返回 (nil, nil)，已包含 FinishRecord 则拒绝续传。
func LoadLocalRebuild(cfg config.Settings, run *runctx.RunContext) (*metrics.RebuildResult, error) {
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
		return nil, fmt.Errorf("%w: run path %q is a directory", ErrInvalidRunFile, runFile)
	}
	res, err := metrics.Rebuild(runFile)
	if err != nil {
		return nil, fmt.Errorf("%w: rebuild local store %q: %w", ErrInvalidRunFile, runFile, err)
	}
	if res.Finished {
		return nil, fmt.Errorf("%w: run %q is already finished, use `swanlab sync` instead", ErrInvalidRunFile, run.RunID)
	}
	return res, nil
}

// ResolveCheckpoint 读取并校验断点文件是否与当前可信本地副本范围一致；
// 缺失、损坏或越界时回退到文件头 (0, HeaderLen, false)。
func ResolveCheckpoint(
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

// CountUploadable 统计 runFile 在 [startOffset, endOffset) 区间内需立即上传的记录条数。
func CountUploadable(runFile string, startOffset, endOffset int64) int64 {
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

// StoreOptions 配置本地 run 文件与异步写盘器打开参数。
type StoreOptions struct {
	SkipStore bool
	Run       *runctx.RunContext
	Resume    *metrics.RebuildResult
	Queue     int
	OnError   func(error)
}

// OpenedStore 描述打开本地存储后的句柄、文件路径与新建标记。
type OpenedStore struct {
	Writer         *store.AsyncWriter
	RunFile        string
	CheckpointFile string
	Created        bool
}

// OpenStore 新建或追加打开本地 run-*.swanlab 并启动异步写盘协程。
func OpenStore(opts StoreOptions) (OpenedStore, error) {
	if opts.SkipStore {
		return OpenedStore{}, nil
	}
	runFile, err := opts.Run.RunFile()
	if err != nil {
		return OpenedStore{}, fmt.Errorf("derive run file: %w", err)
	}
	checkpointFile, err := opts.Run.CheckpointFile()
	if err != nil {
		return OpenedStore{}, fmt.Errorf("derive checkpoint file: %w", err)
	}
	var (
		writer    *store.Writer
		created   bool
		resumeSeq uint64
	)
	if opts.Resume != nil {
		resumeSeq = opts.Resume.Records
		writer, err = store.OpenAppend(runFile, opts.Resume.EndOffset)
		if err != nil {
			return OpenedStore{}, fmt.Errorf("open run file for append: %w", err)
		}
	} else {
		if err = os.MkdirAll(opts.Run.RunDir, 0o755); err != nil {
			return OpenedStore{}, fmt.Errorf("create run dir: %w", err)
		}
		writer, err = store.Create(runFile)
		if err != nil {
			return OpenedStore{}, fmt.Errorf("create run file: %w", err)
		}
		created = true
	}

	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{
		Queue:     opts.Queue,
		ResumeSeq: resumeSeq,
		OnError:   opts.OnError,
	})
	if err != nil {
		_ = writer.Close()
		return OpenedStore{}, fmt.Errorf("start store writer: %w", err)
	}
	return OpenedStore{
		Writer:         async,
		RunFile:        runFile,
		CheckpointFile: checkpointFile,
		Created:        created,
	}, nil
}

// ProgressOptions 配置在线模式的进度跟踪与续传重放起点计算。
type ProgressOptions struct {
	Mode           config.Mode
	RunID          string
	RunFile        string
	CheckpointFile string
	HasStore       bool
	Resume         *metrics.RebuildResult
}

// OnlineProgress 返回在线进度追踪器及续传重放区间起点；非在线或无本地 store 时返回零值结构。
type OnlineProgress struct {
	Tracker           *progress.Tracker
	OnStored          func(*recordv1.Record, uint64)
	ReplayFrom        uint64
	ReplayStartOffset int64
	CheckpointSeq     uint64
	CheckpointValid   bool
}

// InitOnlineProgress 在在线且开启本地 store 时初始化断点水位与待重放范围。
func InitOnlineProgress(opts ProgressOptions) OnlineProgress {
	out := OnlineProgress{ReplayStartOffset: store.HeaderLen}
	if !opts.HasStore || opts.Mode != config.ModeOnline {
		return out
	}
	var baseSeq uint64
	if opts.Resume != nil {
		baseSeq, out.ReplayStartOffset, out.CheckpointValid = ResolveCheckpoint(
			opts.CheckpointFile, opts.RunID, opts.RunFile, opts.Resume,
		)
	}
	out.Tracker = progress.NewTracker(baseSeq)
	out.OnStored = out.Tracker.OnStored
	out.CheckpointSeq = baseSeq
	if opts.Resume != nil && baseSeq < opts.Resume.Records {
		out.ReplayFrom = baseSeq + 1
	}
	return out
}
