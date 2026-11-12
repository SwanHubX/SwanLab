package uploader

import (
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// dispatchOrder 是出批分组的固定顺序：column 必须先于依赖它的数据组送达。
var dispatchOrder = []record.Type{
	record.TypeColumn,
	record.TypeScalar,
	record.TypeMedia,
	record.TypeHistogram,
	record.TypeLog,
	record.TypeSave,
}

// recordGroup 是同一类型的记录集合，records 内保持批内原始顺序。
type recordGroup struct {
	typ     record.Type
	records []*recordv1.Record
}

// recordBuffer 是分发侧独占的缓冲区（不加锁），按 num 原地覆盖并归还旧记录配额。
// total_records 由 intake 在接纳时累加，覆盖时在这里冲减。
type recordBuffer struct {
	budget  *budget.Budget
	tracker *Tracker

	order []int64
	items map[int64]*recordv1.Record
}

// newRecordBuffer 构造空缓冲区。
func newRecordBuffer(b *budget.Budget, tracker *Tracker) *recordBuffer {
	return &recordBuffer{
		budget:  b,
		tracker: tracker,
		items:   make(map[int64]*recordv1.Record),
	}
}

// add 写入一条记录；added 与 replaced 互斥，覆盖时归还旧记录配额并冲减总数。
func (b *recordBuffer) add(rec *recordv1.Record) (added, replaced bool) {
	num := rec.GetNum()
	if old, ok := b.items[num]; ok {
		b.items[num] = rec
		b.budget.Release(1, record.Size(old))
		b.tracker.AdjustTotalRecords(-1)
		return false, true
	}
	b.items[num] = rec
	b.order = append(b.order, num)
	return true, false
}

// len 返回当前条数。
func (b *recordBuffer) len() int { return len(b.order) }

// drain 按插入序返回全部记录并清空缓冲区。
func (b *recordBuffer) drain() []*recordv1.Record {
	if len(b.order) == 0 {
		return nil
	}
	out := make([]*recordv1.Record, 0, len(b.order))
	for _, num := range b.order {
		out = append(out, b.items[num])
	}
	b.order = nil
	b.items = make(map[int64]*recordv1.Record)
	return out
}

// groupRecords 按 dispatchOrder 分组，组内保持批内顺序；未知类型按首次出现顺序追加在末尾。
func groupRecords(records []*recordv1.Record) []recordGroup {
	if len(records) == 0 {
		return nil
	}

	buckets := make(map[record.Type][]*recordv1.Record, len(dispatchOrder))
	seen := make(map[record.Type]struct{}, len(dispatchOrder))
	var extra []record.Type
	for _, rec := range records {
		typ := record.TypeOf(rec)
		buckets[typ] = append(buckets[typ], rec)
		if _, ok := seen[typ]; ok {
			continue
		}
		seen[typ] = struct{}{}
		if !isDispatched(typ) {
			extra = append(extra, typ)
		}
	}

	groups := make([]recordGroup, 0, len(buckets))
	for _, typ := range dispatchOrder {
		if recs := buckets[typ]; len(recs) > 0 {
			groups = append(groups, recordGroup{typ: typ, records: recs})
		}
	}
	for _, typ := range extra {
		groups = append(groups, recordGroup{typ: typ, records: buckets[typ]})
	}
	return groups
}

// isDispatched 报告 typ 是否在固定分发顺序中。
func isDispatched(typ record.Type) bool {
	for _, candidate := range dispatchOrder {
		if candidate == typ {
			return true
		}
	}
	return false
}

// dispatch 按固定类型顺序处理一批；失败时返回该组剩余及全部后续组。
func (u *Uploader) dispatch(records []*recordv1.Record) ([]*recordv1.Record, bool) {
	groups := groupRecords(records)
	var failures []Failure
	for i, group := range groups {
		tail := flattenGroups(groups[i+1:])
		u.pending.Store(int64(len(group.records) + len(tail)))
		sender := u.registry.get(group.typ)
		if sender == nil {
			u.terminateUnsupported(group)
			u.pending.Store(int64(len(tail)))
			continue
		}

		result := safeUpload(u.ctx, sender, group.records)
		if !conserved(group.records, result) {
			result = Result{
				Remaining: group.records,
				Failures: []Failure{{
					Message: "sender result does not account for every record",
					Records: group.records,
				}},
			}
		}
		u.confirm(result.Confirmed)
		failures = append(failures, result.Failures...)

		if len(result.Remaining) > 0 {
			pending := remainingInBatchOrder(group.records, result.Remaining)
			pending = append(pending, tail...)
			u.pending.Store(int64(len(pending)))
			// 随 Confirmed 终结的确定性失败不会再进入下一轮，必须在本轮立即告警；
			// 否则重试成功后该失败会从 LastFailures 消失，用户只能看到通用重试提示。
			u.warnTerminal(terminalFailures(result.Failures, result.Remaining))
			u.setFailures(failures)
			u.warn.Warn()
			return pending, true
		}
		if len(result.Failures) > 0 {
			u.warnTerminal(result.Failures)
		}
		u.pending.Store(int64(len(tail)))
	}
	u.pending.Store(0)
	u.setFailures(failures)
	u.warn.Reset()
	return nil, false
}

// flattenGroups 保持固定组序和组内记录顺序。
func flattenGroups(groups []recordGroup) []*recordv1.Record {
	var records []*recordv1.Record
	for _, group := range groups {
		records = append(records, group.records...)
	}
	return records
}

// confirm 终结记录：归还配额并计入上传进度。
func (u *Uploader) confirm(records []*recordv1.Record) {
	for _, rec := range records {
		u.budget.Release(1, record.Size(rec))
	}
	u.tracker.AdvanceRecords(int64(len(records)))
}

// terminateUnsupported 终结未注册类型的记录：归还配额、计入上传进度并告警一次。
func (u *Uploader) terminateUnsupported(group recordGroup) {
	u.confirm(group.records)
	if _, warned := u.warnedTypes[group.typ]; warned {
		return
	}
	u.warnedTypes[group.typ] = struct{}{}
	u.warnUnsupported(group.typ, len(group.records))
}

// warnTerminal 输出不会被重试的确定性失败，避免静默丢弃数据。
func (u *Uploader) warnTerminal(failures []Failure) {
	for _, failure := range failures {
		u.warnFailure(failure.Message)
	}
}

// setFailures 保存最近一次完整调度尝试的失败集合。
func (u *Uploader) setFailures(failures []Failure) {
	u.failuresMu.Lock()
	u.failures = append(u.failures[:0], failures...)
	u.failuresMu.Unlock()
}
