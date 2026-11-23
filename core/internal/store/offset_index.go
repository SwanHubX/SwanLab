package store

import (
	"sort"
	"sync"
)

// DefaultIndexStride 是 OffsetIndex 的默认锚点间隔（条）。
const DefaultIndexStride = 256

// recentEndCap 是近期记录结束偏移缓存的槽位数（须为 2 的幂，16384 槽约 256 KiB）。
const recentEndCap = 1 << 14

// Anchor 是 storeSeq → 文件偏移的稀疏定位点；Offset 为该条记录头的位置。
type Anchor struct {
	StoreSeq uint64
	Offset   int64
}

// OffsetIndex 按固定间隔记录 (storeSeq, offset)，并缓存近期记录的结束偏移，供断点快速定位。
// 只存在于内存（重启后为空）；登记由 AsyncWriter 单协程串行完成，读取方可并发查询。
type OffsetIndex struct {
	stride  uint64
	mu      sync.RWMutex
	anchors []Anchor
	recent  []Anchor
}

// NewOffsetIndex 创建锚点间隔为 stride 的索引；stride <= 0 时取 DefaultIndexStride。
func NewOffsetIndex(stride int) *OffsetIndex {
	if stride <= 0 {
		stride = DefaultIndexStride
	}
	return &OffsetIndex{stride: uint64(stride)}
}

// Record 缓存第 seq 条记录的结束偏移 end，并在命中锚点间隔时登记起点 start。
// 非单调的稀疏锚点序号会被忽略（AsyncWriter 保证递增）。
func (x *OffsetIndex) Record(seq uint64, start, end int64) {
	if seq == 0 {
		return
	}
	x.mu.Lock()
	defer x.mu.Unlock()
	if x.recent == nil {
		x.recent = make([]Anchor, recentEndCap)
	}
	x.recent[seq&(recentEndCap-1)] = Anchor{StoreSeq: seq, Offset: end}
	x.addAnchorLocked(seq, start)
}

// Add 在 seq 命中锚点间隔时登记 (seq, offset)；非单调序号被忽略（AsyncWriter 保证递增）。
func (x *OffsetIndex) Add(seq uint64, offset int64) {
	if seq == 0 || seq%x.stride != 0 {
		return
	}
	x.mu.Lock()
	defer x.mu.Unlock()
	x.addAnchorLocked(seq, offset)
}

func (x *OffsetIndex) addAnchorLocked(seq uint64, offset int64) {
	if seq%x.stride != 0 {
		return
	}
	if n := len(x.anchors); n > 0 && seq <= x.anchors[n-1].StoreSeq {
		return
	}
	x.anchors = append(x.anchors, Anchor{StoreSeq: seq, Offset: offset})
}

// LookupEnd 在近期环形缓存中查找第 seq 条记录的结束偏移；limit >= 0 时要求结果不超过 limit。
func (x *OffsetIndex) LookupEnd(seq uint64, limit int64) (int64, bool) {
	if x == nil || seq == 0 {
		return 0, false
	}
	x.mu.RLock()
	defer x.mu.RUnlock()
	if len(x.recent) == 0 {
		return 0, false
	}
	slot := x.recent[seq&(recentEndCap-1)]
	if slot.StoreSeq != seq || slot.Offset < HeaderLen {
		return 0, false
	}
	if limit >= 0 && slot.Offset > limit {
		return 0, false
	}
	return slot.Offset, true
}

// LocateFlushed 返回 StoreSeq <= seq 且 Offset < limit 的最后一个锚点，
// 供回填只从已刷盘区域开始读；ok 为 false 时应从文件头开始读。
//
// 锚点由 AsyncWriter 按 (StoreSeq, Offset) 单调登记，因此先二分出落在
// 刷盘上界内的前缀，再在 prefix 内找不晚于 seq 的最后一个锚点。
func (x *OffsetIndex) LocateFlushed(seq uint64, limit int64) (Anchor, bool) {
	x.mu.RLock()
	defer x.mu.RUnlock()
	n := sort.Search(len(x.anchors), func(i int) bool {
		return x.anchors[i].Offset >= limit
	})
	if n == 0 {
		return Anchor{}, false
	}
	i := sort.Search(n, func(i int) bool {
		return x.anchors[i].StoreSeq > seq
	})
	if i == 0 {
		return Anchor{}, false
	}
	return x.anchors[i-1], true
}

// Locate 返回 StoreSeq <= seq 的最后一个锚点；ok 为 false 时应从文件头开始读。
func (x *OffsetIndex) Locate(seq uint64) (Anchor, bool) {
	x.mu.RLock()
	defer x.mu.RUnlock()
	i := sort.Search(len(x.anchors), func(i int) bool {
		return x.anchors[i].StoreSeq > seq
	})
	if i == 0 {
		return Anchor{}, false
	}
	return x.anchors[i-1], true
}

// Len 返回已登记的锚点数。
func (x *OffsetIndex) Len() int {
	x.mu.RLock()
	defer x.mu.RUnlock()
	return len(x.anchors)
}
