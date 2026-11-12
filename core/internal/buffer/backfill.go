// 本文件承载暂停投递期间的本地补投：按调用方给定的顺序把已落盘记录送回上传通道。
package buffer

import (
	"context"
	"errors"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// Replay 把补投记录按给定顺序投入上传通道，并按实际用量占用预算与槽位。
//
// 仅在暂停投递期间可用（存在上传通道的会话）：记录排在实时记录之前，
// 调用方负责按已刷盘边界读取、保证 seq 递增且不重复；本方法不参与接纳计数。
// 预算不足时阻塞到配额释放或 ctx 结束；单条记录超过字节上限时超额占用，
// 与回收补投一致，避免永久阻塞。单批超过上传通道容量返回 ErrBatchTooLarge。
func (t *Intake) Replay(ctx context.Context, records []*recordv1.Record) error {
	if len(records) == 0 {
		return nil
	}
	for _, rec := range records {
		if rec == nil {
			return errors.New("buffer: nil record")
		}
	}

	acquired := 0
	for _, rec := range records {
		size := record.Size(rec)
		err := t.budget.Acquire(ctx, 1, size)
		if errors.Is(err, budget.ErrTooLarge) {
			t.budget.Adjust(1, size)
		} else if err != nil {
			releaseReplayBudget(t.budget, records[:acquired])
			return err
		}
		acquired++
	}

	if err := t.reserveUpload(ctx, len(records)); err != nil {
		releaseReplayBudget(t.budget, records)
		if errors.Is(err, ErrBatchTooLarge) ||
			errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
			return err
		}
		return overload(err)
	}

	t.mu.Lock()
	closed := t.closed
	paused := t.pauseSeq != 0
	if closed || !paused {
		t.releaseUploadLocked(len(records))
		t.mu.Unlock()
		releaseReplayBudget(t.budget, records)
		if closed {
			return ErrClosed
		}
		return ErrNotPaused
	}
	// 槽位已在锁外预留，发送不会阻塞。
	for _, rec := range records {
		t.upload <- rec
	}
	t.releaseUploadLocked(len(records))
	t.mu.Unlock()
	return nil
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

// releaseReplayBudget 归还 Replay 已占用的配额；调用方不得持有 intake 锁。
func releaseReplayBudget(b *budget.Budget, records []*recordv1.Record) {
	for _, rec := range records {
		b.Release(1, record.Size(rec))
	}
}
