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

// PauseState 返回暂停状态与首个未投递记录的 storeSeq（0 表示无）。
func (t *Intake) PauseState() (paused bool, storeSeq uint64) {
	t.mu.Lock()
	defer t.mu.Unlock()
	return t.paused, t.pauseSeq
}

// ResumeUpload 恢复上传投递。调用方必须先补传 pauseSeq 之后的本地记录；
// 恢复后仍可能因旧票据或预算再次进入暂停，须重新观察 PauseState。
func (t *Intake) ResumeUpload() {
	t.mu.Lock()
	t.paused = false
	t.pauseSeq = 0
	t.mu.Unlock()
}

// DrainRecovery 尝试把回收队列投入上传；返回仍未投递的记录数。
// 本地已降级时会把未确认落盘的暂停记录一并转入回收队列，可重复调用。
func (t *Intake) DrainRecovery(ctx context.Context) (int, error) {
	t.mu.Lock()
	if t.degraded.Load() {
		t.collectUnconfirmedLocked()
	}
	t.mu.Unlock()
	for {
		t.mu.Lock()
		if len(t.recovery) == 0 {
			t.mu.Unlock()
			return 0, nil
		}
		if t.closed {
			n := len(t.recovery)
			t.mu.Unlock()
			return n, ErrClosed
		}
		if err := ctx.Err(); err != nil {
			n := len(t.recovery)
			t.mu.Unlock()
			return n, err
		}
		rec := t.recovery[0]
		if !t.sendRecoveryLocked(rec) {
			n := len(t.recovery)
			t.mu.Unlock()
			return n, nil
		}
		t.recovery[0] = nil
		t.recovery = t.recovery[1:]
		if len(t.recovery) == 0 {
			t.recovery = nil
		}
		t.mu.Unlock()
	}
}

// pauseLocked 记录暂停起点、保留未确认记录、暂存 END 保存并归还本批上传预留。
// 若 Resume 之后仍有旧票据提交，这里会开启新的回填区间；需持有 mu。
func (t *Intake) pauseLocked(records []*recordv1.Record, immediate, deferred []int, tk Ticket, seqOf []uint64) {
	if t.pauseSeq == 0 && len(immediate) > 0 {
		t.paused = true
		t.pauseSeq = seqOf[immediate[0]]
	}
	for _, i := range immediate {
		t.unconfirmed = append(t.unconfirmed, unconfirmedRec{seq: seqOf[i], rec: records[i]})
	}
	for _, i := range deferred {
		t.pending = append(t.pending, pendingSave{rec: records[i]})
	}
	if tk.upload > 0 {
		t.budget.Release(int64(tk.upload), tk.bytes)
		t.releaseUploadLocked(tk.upload)
	}
}

// fallbackUploadLocked 处理暂停票据在 Commit 前失去本地副本的情况：
// 有槽位就按条补扣并投递，没有槽位则留在回收队列，避免持锁阻塞。
func (t *Intake) fallbackUploadLocked(records []*recordv1.Record, immediate, deferred []int, tk Ticket) {
	if tk.upload > 0 {
		t.budget.Release(int64(tk.upload), tk.bytes)
		t.releaseUploadLocked(tk.upload)
	}
	for _, i := range immediate {
		if !t.sendRecoveryLocked(records[i]) {
			t.recovery = append(t.recovery, records[i])
		}
	}
	for _, i := range deferred {
		t.pending = append(t.pending, pendingSave{rec: records[i]})
	}
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
