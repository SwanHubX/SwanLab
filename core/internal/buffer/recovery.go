// 本文件承载暂停投递后的未确认记录、回填区间与本地失效后的回收补传。
package buffer

import (
	"context"

	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// unconfirmedRec 是已进入本地写入队列、但尚未确认刷盘的暂停记录。
// 本地写失败时依赖它把记录重新投入上传，避免既无本地副本又不在上传管线。
type unconfirmedRec struct {
	seq uint64
	rec *recordv1.Record
}

// PauseState 返回暂停状态与首个未投递记录的 storeSeq（0 表示未暂停）。
func (t *Intake) PauseState() (paused bool, storeSeq uint64) {
	t.mu.Lock()
	defer t.mu.Unlock()
	return t.pauseSeq != 0, t.pauseSeq
}

// ResumeUpload 恢复上传投递。调用方必须先补传 pauseSeq 之后的本地记录；
// 恢复后仍可能因旧票据或预算再次进入暂停，须重新观察 PauseState。
func (t *Intake) ResumeUpload() {
	t.mu.Lock()
	t.pauseSeq = 0
	t.mu.Unlock()
}

// DrainRecovery 尝试把回收队列投入上传；返回仍未投递的记录数。
// 本地已降级时会把未确认落盘的暂停记录一并转入回收队列，可重复调用。
func (t *Intake) DrainRecovery(ctx context.Context) (int, error) {
	t.mu.Lock()
	defer t.mu.Unlock()
	if t.degraded.Load() {
		t.collectUnconfirmedLocked()
	}
	if len(t.recovery) == 0 {
		return 0, nil
	}
	if t.closed {
		return len(t.recovery), ErrClosed
	}
	if err := ctx.Err(); err != nil {
		return len(t.recovery), err
	}
	t.drainRecoveryLocked()
	return len(t.recovery), nil
}

// drainRecoveryLocked 按 FIFO 尝试补发回收记录；槽位不足时停在最早未投递的
// 记录上，返回是否已全部送出。调用方需持有 mu。
func (t *Intake) drainRecoveryLocked() bool {
	for len(t.recovery) > 0 {
		if !t.sendRecoveryLocked(t.recovery[0]) {
			return false
		}
		t.recovery[0] = nil
		t.recovery = t.recovery[1:]
	}
	t.recovery = nil
	return true
}

// pauseLocked 记录暂停起点、保留未确认记录、暂存 END 保存并归还本批上传槽位。
// 延迟记录从 Reserve 起就占用配额，这里不释放；若 Resume 之后仍有旧票据提交，
// 这里会开启新的回填区间。返回是否新发布了暂停状态。需持有 mu。
func (t *Intake) pauseLocked(records []*recordv1.Record, immediate, deferred []int, tk Ticket, seqOf []uint64) bool {
	published := false
	if t.pauseSeq == 0 && len(immediate) > 0 {
		// 暂停状态只在暂停记录真正提交时发布，票据被放弃不会遗留暂停。
		t.pauseSeq = seqOf[immediate[0]]
		published = true
	}
	for _, i := range immediate {
		t.unconfirmed = append(t.unconfirmed, unconfirmedRec{seq: seqOf[i], rec: records[i]})
	}
	for _, i := range deferred {
		t.pending = append(t.pending, records[i])
	}
	if tk.upload > 0 {
		// 只归还立即记录的配额；延迟记录继续持有到 Close。
		t.budget.Adjust(int64(len(deferred))-int64(tk.upload), deferredBytes(records, deferred)-tk.bytes)
	}
	t.releaseUploadLocked(tk.slots)
	return published
}

// recoveryUploadLocked 处理提交时没有可靠本地副本的批次：并入回收队列，
// 排在已有记录之后按序补传，绝不越过更早的记录。槽位不足时留在队列中，
// 由后续 DrainRecovery 或 Close 继续。需持有 mu。
func (t *Intake) recoveryUploadLocked(records []*recordv1.Record, immediate, deferred []int, tk Ticket) {
	if tk.upload > 0 {
		// 延迟记录继续持有配额；立即记录随回收队列在发送时补扣。
		t.budget.Adjust(int64(len(deferred))-int64(tk.upload), deferredBytes(records, deferred)-tk.bytes)
	}
	for _, i := range immediate {
		t.recovery = append(t.recovery, records[i])
	}
	for _, i := range deferred {
		t.pending = append(t.pending, records[i])
	}
	t.releaseUploadLocked(tk.slots)
	t.drainRecoveryLocked()
}

// deferredBytes 统计延迟记录的序列化字节数。
func deferredBytes(records []*recordv1.Record, deferred []int) int64 {
	var total int64
	for _, i := range deferred {
		total += record.Size(records[i])
	}
	return total
}

// sendRecoveryLocked 为一条回收记录占用槽位、补扣配额并投入上传；需持有 mu。
// 本地副本已失效，因此不再等待预算，超限部分由回填积压量约束。
func (t *Intake) sendRecoveryLocked(rec *recordv1.Record) bool {
	if !t.tryReserveUploadLocked(1) {
		return false
	}
	t.budget.Adjust(1, record.Size(rec))
	t.upload <- rec
	t.releaseUploadLocked(1)
	return true
}

// trimUnconfirmedLocked 丢弃已确认刷盘的前缀；调用方需持有 mu。
func (t *Intake) trimUnconfirmedLocked() {
	if len(t.unconfirmed) == 0 {
		return
	}
	flushed := t.flushedSeq()
	i := 0
	for i < len(t.unconfirmed) && t.unconfirmed[i].seq <= flushed {
		t.unconfirmed[i] = unconfirmedRec{}
		i++
	}
	if i > 0 {
		t.unconfirmed = t.unconfirmed[i:]
	}
	if len(t.unconfirmed) == 0 {
		t.unconfirmed = nil
	}
}

// collectUnconfirmedLocked 把未确认落盘的暂停记录转入回收队列；调用方需持有 mu。
func (t *Intake) collectUnconfirmedLocked() {
	if len(t.unconfirmed) == 0 {
		return
	}
	flushed := t.flushedSeq()
	for _, u := range t.unconfirmed {
		if u.seq > flushed {
			t.recovery = append(t.recovery, u.rec)
		}
	}
	t.unconfirmed = nil
}

func (t *Intake) flushedSeq() uint64 {
	if t.flushed == nil {
		return 0
	}
	return t.flushed()
}

// restoreRecovery 把未发出的回收记录放回队头，保持原有顺序。
func (t *Intake) restoreRecovery(items []*recordv1.Record) {
	if len(items) == 0 {
		return
	}
	t.mu.Lock()
	rest := append([]*recordv1.Record(nil), items...)
	t.recovery = append(rest, t.recovery...)
	t.mu.Unlock()
}
