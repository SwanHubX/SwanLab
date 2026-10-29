package uploader

import (
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// dispatchOrder 是出批分组的固定顺序；门控由轨道依赖表决定，与此无关。
var dispatchOrder = []record.Type{
	record.TypeColumn,
	record.TypeScalar,
	record.TypeMedia,
	record.TypeLog,
	record.TypeSave,
}

// recordGroup 是同一类型的记录集合，records 内保持批内原始顺序。
type recordGroup struct {
	typ     record.Type
	records []*recordv1.Record
}

// recordBuffer 是 batcher 独占的缓冲区（不加锁），按 num 原地覆盖并归还旧记录配额。
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

// add 写入一条记录；added 与 replaced 互斥，覆盖时归还旧记录配额。
func (b *recordBuffer) add(rec *recordv1.Record) (added, replaced bool) {
	num := rec.GetNum()
	if old, ok := b.items[num]; ok {
		b.items[num] = rec
		b.budget.Release(1, record.Size(old))
		return false, true
	}
	b.items[num] = rec
	b.order = append(b.order, num)
	b.tracker.AddTotalRecords(1)
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
