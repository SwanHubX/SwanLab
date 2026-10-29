package uploader

import (
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// dispatchOrder 是出批分组的固定顺序：列记录永远排在最前，其余类型顺序固定，
// 保证日志与测试可读。跨轨门控由轨道依赖表决定，不依赖这里的顺序。
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

// recordBuffer 是攒批缓冲区，由 batcher 协程独占，因此不加锁。
//
// 以 num 为主键原地覆盖：同一 num 再次写入时保留首次出现的位置，只替换
// 记录内容，并立刻归还旧记录占用的配额；新位置才计入总记录数。
type recordBuffer struct {
	budget  *budget.Budget
	tracker *Tracker

	order []int64
	items map[int64]*recordv1.Record
}

// newRecordBuffer 构造空缓冲区；budget 与 tracker 必须非空。
func newRecordBuffer(b *budget.Budget, tracker *Tracker) *recordBuffer {
	return &recordBuffer{
		budget:  b,
		tracker: tracker,
		items:   make(map[int64]*recordv1.Record),
	}
}

// add 写入一条记录。added 表示新增了一个位置，replaced 表示覆盖了同 num 的
// 旧记录；两者互斥。覆盖时旧记录的配额在这里归还。
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

// len 返回当前条数，覆盖不改变长度。
func (b *recordBuffer) len() int { return len(b.order) }

// drain 按插入序返回全部记录并清空缓冲区；items 换新 map，避免长期持有大对象。
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

// groupRecords 按 dispatchOrder 分组，组内保持批内原始顺序。
// dispatchOrder 之外的类型（未知 oneof、start/finish）按首次出现顺序追加在
// 末尾，由调用方决定如何终结。
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

// isDispatched 报告 typ 是否属于固定分发顺序中的上传类型。
func isDispatched(typ record.Type) bool {
	for _, candidate := range dispatchOrder {
		if candidate == typ {
			return true
		}
	}
	return false
}
