// 本文件承载暂停投递期间的本地补投：按调用方给定的顺序把已落盘记录送回上传通道。
package buffer

import (
	"context"
	"errors"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// replayChunkMax 是单次批量投递的最大记录数，避免长时间持有 intake 锁。
const replayChunkMax = 256

// Replay 按顺序把补投记录投入上传通道，返回成功投递的条数。
//
// 仅在暂停投递期间可用（存在上传通道的会话）：记录排在实时记录之前，
// 调用方负责按已刷盘边界读取、保证 seq 递增且不重复；本方法不参与接纳计数。
// 预算按当前可用容量批量占用，不足时等待到配额释放或 ctx 结束；ctx 结束或
// 投递中断时返回已投递条数，调用方按返回值推进游标，未投递记录可安全重试。
func (t *Intake) Replay(ctx context.Context, records []*recordv1.Record) (int, error) {
	for _, rec := range records {
		if rec == nil {
			return 0, errors.New("buffer: nil record")
		}
	}

	sent := 0
	for sent < len(records) {
		batch := records[sent:]
		n := t.replayChunk(batch)
		bytes := replayBytes(batch[:n])
		for n > 1 && !t.budget.TryAcquire(int64(n), bytes) {
			n /= 2
			bytes = replayBytes(batch[:n])
		}
		if n == 1 {
			// 单条记录超过字节上限：与回收补投一致，超额占用避免永久阻塞。
			err := t.budget.Acquire(ctx, 1, bytes)
			if errors.Is(err, budget.ErrTooLarge) {
				t.budget.Adjust(1, bytes)
			} else if err != nil {
				return sent, err
			}
		}

		if err := t.reserveUpload(ctx, n); err != nil {
			t.budget.Release(int64(n), bytes)
			if errors.Is(err, ErrBatchTooLarge) ||
				errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
				return sent, err
			}
			return sent, overload(err)
		}

		t.mu.Lock()
		closed := t.closed
		paused := t.pauseSeq != 0
		if !closed && paused {
			for _, rec := range batch[:n] {
				t.upload <- rec
			}
			sent += n
		}
		t.releaseUploadLocked(n)
		t.mu.Unlock()
		if closed || !paused {
			t.budget.Release(int64(n), bytes)
			if closed {
				return sent, ErrClosed
			}
			return sent, ErrNotPaused
		}
	}
	return sent, nil
}

// replayChunk 返回单批可尝试占用的记录数，不超过预算条数上限。
func (t *Intake) replayChunk(records []*recordv1.Record) int {
	n := min(len(records), replayChunkMax)
	if limit := t.budget.Limits().Records; limit > 0 && int64(n) > limit {
		n = int(limit)
	}
	return n
}

// replayBytes 统计一批记录的序列化字节数。
func replayBytes(records []*recordv1.Record) int64 {
	var total int64
	for _, rec := range records {
		total += record.Size(rec)
	}
	return total
}

// Handover 在锁内补齐补投尾部并原子恢复实时投递。
//
// sync 先把已提交记录全部刷盘（nil 表示调用方已确认刷盘完成）；fill 随后
// 返回尾部记录与目标 storeSeq。记录必须按序、已落盘且不重复，目标必须等于
// 当前最后提交的 storeSeq。失败时不发送任何记录并保留暂停状态。
//
// 成功返回后 pauseSeq 已清零，其后提交的记录直接进入上传通道。
func (t *Intake) Handover(
	ctx context.Context,
	sync func(context.Context) error,
	fill func() ([]*recordv1.Record, uint64, error),
) error {
	if fill == nil {
		return errors.New("buffer: nil backfill filler")
	}

	t.mu.Lock()
	defer t.mu.Unlock()
	if t.closed {
		return ErrClosed
	}
	if t.pauseSeq == 0 {
		return ErrNotPaused
	}
	if sync != nil {
		if err := sync(ctx); err != nil {
			return err
		}
	}

	records, target, err := fill()
	if err != nil {
		return err
	}
	if target != t.nextSeq {
		return ErrBackfillStale
	}
	for _, rec := range records {
		if rec == nil {
			return errors.New("buffer: nil record")
		}
	}
	if len(records) == 0 {
		t.pauseSeq = 0
		return nil
	}
	if !t.tryReserveUploadLocked(len(records)) {
		return ErrBackfillBlocked
	}
	for _, rec := range records {
		// 暂停期间这些记录没有占用预算；尾部长度由调用方限定，超额部分有界。
		t.budget.Adjust(1, record.Size(rec))
	}
	// 槽位已在锁内预留，发送不会阻塞。
	for _, rec := range records {
		t.upload <- rec
	}
	t.releaseUploadLocked(len(records))
	t.pauseSeq = 0
	return nil
}
