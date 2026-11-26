// Package progress 追踪已落盘记录的连续上传确认水位（store_seq）。
package progress

import (
	"sync"

	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// Tracker 追踪已落盘记录的连续上传确认水位（store_seq）。
//
// 所有方法并发安全，分别由 intake.OnStored、backfill 扫描与 uploader.OnConfirm / OnDrop 调用。
// Start / Finish / SAVE_POLICY_END+CUSTOM 等实时阶段不经过 uploader 的记录会在注册时
// 立即视作已完结，以免阻塞连续水位前进。
type Tracker struct {
	mu       sync.Mutex
	next     int64                        // 下一个等待完结的 storeSeq；连续确认水位为 next - 1
	ptrs     map[*recordv1.Record]int64   // 活跃记录指针 -> storeSeq
	seqPtrs  map[int64][]*recordv1.Record // storeSeq -> 活跃指针列表，保留同序号多指针直到任一确认
	resolved map[int64]struct{}           // 已完结但前面仍有缺口的 storeSeq
}

// NewTracker 构造以 baseSeq 为初始连续水位的进度追踪器。
func NewTracker(baseSeq int64) *Tracker {
	return &Tracker{
		next:     baseSeq + 1,
		ptrs:     make(map[*recordv1.Record]int64),
		seqPtrs:  make(map[int64][]*recordv1.Record),
		resolved: make(map[int64]struct{}),
	}
}

// OnStored 注册一条已分配 storeSeq 的记录。
// 同一序号重复注册不同指针时保留旧指针映射，任一指针确认即清理同一序号关联的所有指针。
// 非实时上传类型会立即视作已完结。
func (p *Tracker) OnStored(rec *recordv1.Record, seq int64) {
	if p == nil || rec == nil || seq <= 0 {
		return
	}
	p.mu.Lock()
	defer p.mu.Unlock()

	if seq < p.next {
		return
	}
	if !record.IsImmediateUploadable(rec) {
		p.clearSeqPtrsLocked(seq)
		p.resolveSeqLocked(seq)
		return
	}
	if _, done := p.resolved[seq]; done {
		return
	}
	if existing, ok := p.ptrs[rec]; ok && existing == seq {
		return
	}
	p.ptrs[rec] = seq
	p.seqPtrs[seq] = append(p.seqPtrs[seq], rec)
}

// Resolve 标记一批记录已完成上传确认或确定性终结。
// 未注册过的记录（例如崩溃 error_log 等仅上传不落盘记录）会被忽略。
func (p *Tracker) Resolve(records []*recordv1.Record) {
	if p == nil || len(records) == 0 {
		return
	}
	p.mu.Lock()
	defer p.mu.Unlock()
	for _, rec := range records {
		p.resolvePtrLocked(rec)
	}
}

// ResolveOne 标记单条被丢弃或终结的记录已完成（例如同 num 覆盖）。
func (p *Tracker) ResolveOne(rec *recordv1.Record) {
	if p == nil || rec == nil {
		return
	}
	p.mu.Lock()
	defer p.mu.Unlock()
	p.resolvePtrLocked(rec)
}

// Watermark 返回已连续完结的最大 storeSeq。
func (p *Tracker) Watermark() int64 {
	if p == nil {
		return 0
	}
	p.mu.Lock()
	w := p.next - 1
	p.mu.Unlock()
	return w
}

func (p *Tracker) resolvePtrLocked(rec *recordv1.Record) {
	if rec == nil {
		return
	}
	seq, ok := p.ptrs[rec]
	if !ok {
		return
	}
	p.clearSeqPtrsLocked(seq)
	p.resolveSeqLocked(seq)
}

func (p *Tracker) clearSeqPtrsLocked(seq int64) {
	for _, ptr := range p.seqPtrs[seq] {
		delete(p.ptrs, ptr)
	}
	delete(p.seqPtrs, seq)
}

func (p *Tracker) resolveSeqLocked(seq int64) {
	if seq < p.next {
		return
	}
	p.resolved[seq] = struct{}{}
	for {
		if _, ok := p.resolved[p.next]; !ok {
			return
		}
		delete(p.resolved, p.next)
		p.next++
	}
}
