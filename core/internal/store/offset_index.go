package store

import (
	"sort"
	"sync"
)

// DefaultIndexStride 是 OffsetIndex 的默认锚点间隔（条）。
const DefaultIndexStride = 256

// Anchor 是 storeSeq → 文件偏移的稀疏定位点；Offset 为该条记录头的位置。
type Anchor struct {
	StoreSeq uint64
	Offset   int64
}

// OffsetIndex 按固定间隔记录 (storeSeq, offset)，让回填从任意序号附近开始读取。
// 只存在于内存（重启后为空）；登记由 AsyncWriter 单协程串行完成，读取方可并发 Locate。
type OffsetIndex struct {
	stride  uint64
	mu      sync.RWMutex
	anchors []Anchor
}

// NewOffsetIndex 创建锚点间隔为 stride 的索引；stride <= 0 时取 DefaultIndexStride。
func NewOffsetIndex(stride int) *OffsetIndex {
	if stride <= 0 {
		stride = DefaultIndexStride
	}
	return &OffsetIndex{stride: uint64(stride)}
}

// Add 在 seq 命中锚点间隔时登记 (seq, offset)；非单调序号被忽略（AsyncWriter 保证递增）。
func (x *OffsetIndex) Add(seq uint64, offset int64) {
	if seq == 0 || seq%x.stride != 0 {
		return
	}
	x.mu.Lock()
	defer x.mu.Unlock()
	if n := len(x.anchors); n > 0 && seq <= x.anchors[n-1].StoreSeq {
		return
	}
	x.anchors = append(x.anchors, Anchor{StoreSeq: seq, Offset: offset})
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
