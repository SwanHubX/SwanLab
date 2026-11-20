package store

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"time"
)

// CheckpointVersion 是 checkpoint 文件格式版本号。
const CheckpointVersion = 1

// Checkpoint 记录已上传确认前缀及其在 run-*.swanlab 中的结束偏移。
type Checkpoint struct {
	Version   int       `json:"version"`
	RunID     string    `json:"run_id"`
	StoreSeq  uint64    `json:"store_seq"`
	Offset    int64     `json:"offset"`
	UpdatedAt time.Time `json:"updated_at"`
}

// LoadCheckpoint 读取并校验 checkpoint 文件。
// 文件不存在返回 (zero, false, nil)；版本不匹配、JSON 损坏、run_id 不符或偏移非法时
// 返回 (zero, false, err)，调用方按无 checkpoint 处理并回退到从头重放。
func LoadCheckpoint(path, runID string) (Checkpoint, bool, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return Checkpoint{}, false, nil
		}
		return Checkpoint{}, false, err
	}
	var cp Checkpoint
	if err := json.Unmarshal(data, &cp); err != nil {
		return Checkpoint{}, false, fmt.Errorf("store: invalid checkpoint json: %w", err)
	}
	if cp.Version != CheckpointVersion {
		return Checkpoint{}, false, fmt.Errorf("store: checkpoint version %d != %d", cp.Version, CheckpointVersion)
	}
	if cp.RunID == "" {
		return Checkpoint{}, false, errors.New("store: checkpoint missing run_id")
	}
	if runID != "" && cp.RunID != runID {
		return Checkpoint{}, false, fmt.Errorf("store: checkpoint run_id %q != %q", cp.RunID, runID)
	}
	if cp.Offset < HeaderLen {
		return Checkpoint{}, false, fmt.Errorf("store: checkpoint offset %d < header %d", cp.Offset, HeaderLen)
	}
	if cp.StoreSeq > 0 && cp.Offset == HeaderLen {
		return Checkpoint{}, false, fmt.Errorf("store: checkpoint seq %d has header-only offset %d", cp.StoreSeq, cp.Offset)
	}
	return cp, true, nil
}

// WriteCheckpoint 将 checkpoint 原子写入 path：同目录临时文件 -> fsync -> rename -> fsync 目录。
// 任一步失败均不会破坏已有目标文件，并清理临时文件。
func WriteCheckpoint(path string, cp Checkpoint) error {
	return writeCheckpoint(path, cp, writeTempFile)
}

func writeCheckpoint(path string, cp Checkpoint, writeTemp func(*os.File, []byte) error) error {
	if cp.Version == 0 {
		cp.Version = CheckpointVersion
	} else if cp.Version != CheckpointVersion {
		return fmt.Errorf("store: unsupported checkpoint version %d", cp.Version)
	}
	if cp.RunID == "" {
		return errors.New("store: checkpoint missing run_id")
	}
	if cp.Offset < HeaderLen {
		return fmt.Errorf("store: checkpoint offset %d < header %d", cp.Offset, HeaderLen)
	}
	if cp.UpdatedAt.IsZero() {
		cp.UpdatedAt = time.Now().UTC()
	} else {
		cp.UpdatedAt = cp.UpdatedAt.UTC()
	}

	data, err := json.Marshal(cp)
	if err != nil {
		return fmt.Errorf("store: marshal checkpoint: %w", err)
	}

	dir := filepath.Dir(path)
	tmp, err := os.CreateTemp(dir, "."+filepath.Base(path)+".tmp.*")
	if err != nil {
		return fmt.Errorf("store: create temp checkpoint in %s: %w", dir, err)
	}
	tmpName := tmp.Name()
	renamed := false
	defer func() {
		if !renamed {
			_ = os.Remove(tmpName)
		}
	}()

	if err = writeTemp(tmp, data); err != nil {
		_ = tmp.Close()
		return err
	}
	if err = tmp.Sync(); err != nil {
		_ = tmp.Close()
		return fmt.Errorf("store: sync temp checkpoint: %w", err)
	}
	if err = tmp.Close(); err != nil {
		return fmt.Errorf("store: close temp checkpoint: %w", err)
	}
	if err = os.Chmod(tmpName, 0o644); err != nil {
		return fmt.Errorf("store: chmod temp checkpoint: %w", err)
	}
	if err = os.Rename(tmpName, path); err != nil {
		return fmt.Errorf("store: commit checkpoint: %w", err)
	}
	renamed = true
	syncDir(dir)
	return nil
}

func writeTempFile(fp *os.File, data []byte) error {
	if _, err := fp.Write(data); err != nil {
		return fmt.Errorf("store: write temp checkpoint: %w", err)
	}
	return nil
}

// syncDir 持久化目录项，使 rename 结果落盘；失败不影响写入结果。
func syncDir(dir string) {
	d, err := os.Open(dir)
	if err != nil {
		return
	}
	defer func() { _ = d.Close() }()
	_ = d.Sync()
}

// EndOffset 定位第 seq 条记录之后的物理文件偏移（下一条记录的起点）。
// 若 idx 非 nil，先查找不晚于 seq 的稀疏锚点（limit >= 0 时同时约束在 limit 之内），
// 再从锚点起线性扫描至多 stride 条记录；无锚点时从文件头起扫。
func EndOffset(path string, idx *OffsetIndex, seq uint64, limit int64) (int64, error) {
	if seq == 0 {
		return HeaderLen, nil
	}
	var cur uint64
	start := int64(HeaderLen)
	if idx != nil {
		var (
			anchor Anchor
			ok     bool
		)
		if limit >= 0 {
			anchor, ok = idx.LocateFlushed(seq, limit)
		} else {
			anchor, ok = idx.Locate(seq)
		}
		if ok && anchor.StoreSeq > 0 {
			cur = anchor.StoreSeq - 1
			start = anchor.Offset
		}
	}

	r, err := OpenOffset(path, start, limit)
	if err != nil {
		return 0, err
	}
	defer func() { _ = r.Close() }()

	for cur < seq {
		if _, scanErr := r.Scan(); scanErr != nil {
			if errors.Is(scanErr, io.EOF) {
				return 0, fmt.Errorf("store: %w: reached EOF at seq %d before %d", ErrTruncated, cur, seq)
			}
			return 0, scanErr
		}
		cur++
	}
	return r.Position(), nil
}
