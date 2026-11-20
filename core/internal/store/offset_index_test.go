package store

import (
	"bytes"
	"sync"
	"testing"
)

func TestOffsetIndexLocate(t *testing.T) {
	idx := NewOffsetIndex(3)
	for seq := uint64(1); seq <= 9; seq++ {
		idx.Add(seq, int64(seq)*100)
	}
	if got := idx.Len(); got != 3 {
		t.Fatalf("anchors = %d, want 3", got)
	}

	cases := []struct {
		seq  uint64
		want Anchor
		ok   bool
	}{
		{0, Anchor{}, false},
		{1, Anchor{}, false},
		{3, Anchor{StoreSeq: 3, Offset: 300}, true},
		{5, Anchor{StoreSeq: 3, Offset: 300}, true},
		{6, Anchor{StoreSeq: 6, Offset: 600}, true},
		{9, Anchor{StoreSeq: 9, Offset: 900}, true},
		{100, Anchor{StoreSeq: 9, Offset: 900}, true},
	}
	for _, c := range cases {
		got, ok := idx.Locate(c.seq)
		if ok != c.ok || got != c.want {
			t.Fatalf("Locate(%d) = (%+v, %v), want (%+v, %v)", c.seq, got, ok, c.want, c.ok)
		}
	}
}

// LocateFlushed 只能返回已刷盘区域内的锚点：Offset 恰好等于上界的锚点、
// 以及上界之前没有锚点时都必须退回文件头读取。
func TestOffsetIndexLocateFlushed(t *testing.T) {
	idx := NewOffsetIndex(2)
	for seq := uint64(1); seq <= 8; seq++ {
		idx.Add(seq, int64(seq)*100) // 锚点位于 2 / 4 / 6 / 8
	}

	cases := []struct {
		seq   uint64
		limit int64
		want  Anchor
		ok    bool
	}{
		{8, 801, Anchor{StoreSeq: 8, Offset: 800}, true},
		{8, 800, Anchor{StoreSeq: 6, Offset: 600}, true},
		{5, 700, Anchor{StoreSeq: 4, Offset: 400}, true},
		{7, 700, Anchor{StoreSeq: 6, Offset: 600}, true},
		{5, 300, Anchor{StoreSeq: 2, Offset: 200}, true},
		{5, 200, Anchor{}, false},
		{1, 900, Anchor{}, false},
		{0, 900, Anchor{}, false},
		{8, HeaderLen, Anchor{}, false},
	}
	for _, c := range cases {
		got, ok := idx.LocateFlushed(c.seq, c.limit)
		if ok != c.ok || got != c.want {
			t.Fatalf("LocateFlushed(%d, %d) = (%+v, %v), want (%+v, %v)",
				c.seq, c.limit, got, ok, c.want, c.ok)
		}
	}
}

// 未刷盘时不能从任何锚点开始读；Flush 之后锚点才进入可回填区间。
func TestAsyncWriterLocateFlushed(t *testing.T) {
	a, _ := newTestAsync(t, AsyncOptions{Queue: 32, IndexStride: 2})
	for i := 1; i <= 6; i++ {
		a.In() <- Request{StoreSeq: uint64(i), Payload: anyPayload(i)}
	}

	flushedSeq, flushedOffset := a.Flushed()
	if flushedSeq != 0 || flushedOffset != HeaderLen {
		t.Fatalf("Flushed before sync = (%d, %d), want (0, %d)", flushedSeq, flushedOffset, HeaderLen)
	}
	if _, ok := a.Index().LocateFlushed(6, flushedOffset); ok {
		t.Fatal("LocateFlushed before flush should not expose an anchor")
	}

	if err := a.FlushForTest(); err != nil {
		t.Fatalf("FlushForTest: %v", err)
	}
	flushedSeq, flushedOffset = a.Flushed()
	if flushedSeq != 6 {
		t.Fatalf("Flushed seq after sync = %d, want 6", flushedSeq)
	}
	anchor, ok := a.Index().LocateFlushed(flushedSeq, flushedOffset)
	if !ok || anchor.StoreSeq != 6 {
		t.Fatalf("LocateFlushed(%d, %d) = (%+v, %v), want anchor 6", flushedSeq, flushedOffset, anchor, ok)
	}
	if _, ok := a.Index().LocateFlushed(1, flushedOffset); ok {
		t.Fatal("seq before the first anchor must fall back to the file header")
	}
	if err := a.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

func TestOffsetIndexDefaultStride(t *testing.T) {
	idx := NewOffsetIndex(0)
	idx.Add(DefaultIndexStride-1, 1)
	idx.Add(DefaultIndexStride, 2)
	if got := idx.Len(); got != 1 {
		t.Fatalf("anchors = %d, want 1", got)
	}
	got, ok := idx.Locate(DefaultIndexStride)
	if !ok || got.Offset != 2 {
		t.Fatalf("Locate(%d) = (%+v, %v), want offset 2", DefaultIndexStride, got, ok)
	}
}

func TestOffsetIndexIgnoresNonMonotonic(t *testing.T) {
	idx := NewOffsetIndex(1)
	idx.Add(5, 50)
	idx.Add(3, 30) // 乱序
	idx.Add(5, 55) // 重复
	idx.Add(6, 60)

	if got := idx.Len(); got != 2 {
		t.Fatalf("anchors = %d, want 2", got)
	}
	if got, _ := idx.Locate(5); got.Offset != 50 {
		t.Fatalf("Locate(5).Offset = %d, want 50", got.Offset)
	}
}

func TestOffsetIndexConcurrent(t *testing.T) {
	idx := NewOffsetIndex(1)
	const n = 2000

	var wg sync.WaitGroup
	wg.Add(1)
	go func() {
		defer wg.Done()
		for seq := uint64(1); seq <= n; seq++ {
			idx.Add(seq, int64(seq))
		}
	}()
	for r := range 4 {
		wg.Add(1)
		go func(r int) {
			defer wg.Done()
			for i := range n {
				seq := uint64((i*7 + r) % n)
				if got, ok := idx.Locate(seq); ok && got.StoreSeq > seq {
					t.Errorf("Locate(%d) = %d, ahead of requested seq", seq, got.StoreSeq)
					return
				}
				_ = idx.Len()
			}
		}(r)
	}
	wg.Wait()

	if got := idx.Len(); got != n {
		t.Fatalf("anchors = %d, want %d", got, n)
	}
	if got, ok := idx.Locate(n); !ok || got.StoreSeq != n {
		t.Fatalf("Locate(%d) = (%+v, %v)", n, got, ok)
	}
}

// AsyncWriter 必须登记记录头的真实偏移，供 OpenOffset 直接定位。
func TestAsyncWriterIndexAnchors(t *testing.T) {
	a, path := newTestAsync(t, AsyncOptions{Queue: 64, IndexStride: 3})

	const n = 10
	for i := range n {
		a.In() <- Request{StoreSeq: uint64(i + 1), Payload: anyPayload(i)}
	}
	if err := a.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	idx := a.Index()
	if got := idx.Len(); got != 3 {
		t.Fatalf("anchors = %d, want 3", got)
	}
	if _, ok := idx.Locate(2); ok {
		t.Fatal("Locate(2) should fall before the first anchor")
	}
	for _, seq := range []uint64{3, 6, 9} {
		anchor, ok := idx.Locate(seq)
		if !ok || anchor.StoreSeq != seq {
			t.Fatalf("Locate(%d) = (%+v, %v)", seq, anchor, ok)
		}
		got := readLimit(t, path, anchor.Offset, -1)
		if want := n - int(seq) + 1; len(got) != want {
			t.Fatalf("seq %d: read %d records from anchor, want %d", seq, len(got), want)
		}
		if !bytes.Equal(got[0], anyPayload(int(seq)-1)) {
			t.Fatalf("seq %d: anchor points at the wrong record", seq)
		}
	}
}

// 上一条记录写完后块内只剩 3 字节时，锚点必须是补零之后的记录头，而不是填充起点。
func TestAsyncWriterIndexAnchorAfterPadding(t *testing.T) {
	a, path := newTestAsync(t, AsyncOptions{Queue: 8, IndexStride: 1})

	a.In() <- Request{StoreSeq: 1, Payload: payloadN(BlockLen - 2*HeaderLen - 3)}
	a.In() <- Request{StoreSeq: 2, Payload: payloadN(16)}
	if err := a.FlushForTest(); err != nil {
		t.Fatalf("FlushForTest: %v", err)
	}

	anchor, ok := a.Index().Locate(2)
	if !ok {
		t.Fatal("Locate(2) missing")
	}
	if anchor.Offset != BlockLen {
		t.Fatalf("anchor offset = %d, want record head %d", anchor.Offset, int64(BlockLen))
	}
	got := readLimit(t, path, anchor.Offset, -1)
	if len(got) != 1 || !bytes.Equal(got[0], payloadN(16)) {
		t.Fatalf("read from anchor = %d records, want the second record only", len(got))
	}
	if err := a.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}
