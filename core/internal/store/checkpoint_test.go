package store

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"testing"
	"time"
)

func TestCheckpointRoundTrip(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run-r1.swanlab.checkpoint.json")
	ts := time.Date(2026, 11, 13, 2, 3, 4, 123456789, time.UTC)
	want := Checkpoint{
		Version:   CheckpointVersion,
		RunID:     "run-1",
		StoreSeq:  42,
		Offset:    4096,
		UpdatedAt: ts,
	}
	if err := WriteCheckpoint(path, want); err != nil {
		t.Fatalf("WriteCheckpoint: %v", err)
	}

	got, ok, err := LoadCheckpoint(path, "run-1")
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint = (%+v, %v, %v), want ok=true, err=nil", got, ok, err)
	}
	if got.Version != want.Version || got.RunID != want.RunID ||
		got.StoreSeq != want.StoreSeq || got.Offset != want.Offset ||
		!got.UpdatedAt.Equal(want.UpdatedAt) {
		t.Fatalf("LoadCheckpoint = %+v, want %+v", got, want)
	}
}

func TestCheckpointMissingAndInvalid(t *testing.T) {
	dir := t.TempDir()
	missing := filepath.Join(dir, "missing.checkpoint.json")
	if _, ok, err := LoadCheckpoint(missing, "run-1"); err != nil || ok {
		t.Fatalf("missing checkpoint = (ok=%v, err=%v), want (false, nil)", ok, err)
	}

	cases := []struct {
		name string
		raw  string
	}{
		{"corrupt json", `{broken`},
		{"wrong version", `{"version":2,"run_id":"run-1","store_seq":1,"offset":64}`},
		{"missing run_id", `{"version":1,"run_id":"","store_seq":1,"offset":64}`},
		{"wrong run_id", `{"version":1,"run_id":"other","store_seq":1,"offset":64}`},
		{"offset before header", `{"version":1,"run_id":"run-1","store_seq":0,"offset":3}`},
		{"negative store_seq", `{"version":1,"run_id":"run-1","store_seq":-1,"offset":64}`},
		{"nonzero seq at header", `{"version":1,"run_id":"run-1","store_seq":1,"offset":7}`},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			p := filepath.Join(dir, "cp.json")
			if err := os.WriteFile(p, []byte(c.raw), 0o644); err != nil {
				t.Fatalf("WriteFile: %v", err)
			}
			if _, ok, err := LoadCheckpoint(p, "run-1"); ok || err == nil {
				t.Fatalf("LoadCheckpoint = (ok=%v, err=%v), want (false, non-nil)", ok, err)
			}
		})
	}
}

func TestCheckpointAtomicReplacePreservesOldOnFailure(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "run.checkpoint.json")
	oldCP := Checkpoint{
		Version:   CheckpointVersion,
		RunID:     "run-1",
		StoreSeq:  10,
		Offset:    512,
		UpdatedAt: time.Unix(1700000000, 0).UTC(),
	}
	if err := WriteCheckpoint(path, oldCP); err != nil {
		t.Fatalf("WriteCheckpoint initial: %v", err)
	}
	oldRaw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile initial: %v", err)
	}

	// 模拟写入临时文件失败：旧文件必须保持原样且不残留临时文件。
	newCP := Checkpoint{
		Version:   CheckpointVersion,
		RunID:     "run-1",
		StoreSeq:  20,
		Offset:    1024,
		UpdatedAt: time.Unix(1700000100, 0).UTC(),
	}
	failWrite := func(fp *os.File, data []byte) error {
		_, _ = fp.Write(data[:len(data)/2])
		return errors.New("injected disk full")
	}
	if werr := writeCheckpoint(path, newCP, true, failWrite); werr == nil {
		t.Fatal("writeCheckpoint with failing temp write should fail")
	}

	afterRaw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile after failure: %v", err)
	}
	if string(afterRaw) != string(oldRaw) {
		t.Fatalf("old checkpoint modified after failure:\ngot:  %s\nwant: %s", afterRaw, oldRaw)
	}
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("ReadDir: %v", err)
	}
	if len(entries) != 1 || entries[0].Name() != filepath.Base(path) {
		t.Fatalf("unexpected temp files left in %s: %v", dir, entries)
	}

	// 覆盖写入成功场景（含运行期 WriteCheckpointFast 与完整同步 WriteCheckpoint）。
	fastCP := newCP
	fastCP.StoreSeq = 15
	fastCP.Offset = 768
	if werr := WriteCheckpointFast(path, fastCP); werr != nil {
		t.Fatalf("WriteCheckpointFast replace: %v", werr)
	}
	gotFast, ok, err := LoadCheckpoint(path, "run-1")
	if err != nil || !ok || gotFast.StoreSeq != fastCP.StoreSeq || gotFast.Offset != fastCP.Offset {
		t.Fatalf("fast checkpoint = (%+v, %v, %v), want seq=%d offset=%d",
			gotFast, ok, err, fastCP.StoreSeq, fastCP.Offset)
	}

	if werr := WriteCheckpoint(path, newCP); werr != nil {
		t.Fatalf("WriteCheckpoint replace: %v", werr)
	}
	got, ok, err := LoadCheckpoint(path, "run-1")
	if err != nil || !ok || got.StoreSeq != newCP.StoreSeq || got.Offset != newCP.Offset {
		t.Fatalf("replaced checkpoint = (%+v, %v, %v), want seq=%d offset=%d",
			got, ok, err, newCP.StoreSeq, newCP.Offset)
	}

	// 非法输入也不能破坏旧文件。
	if werr := WriteCheckpoint(path, Checkpoint{RunID: "", Offset: HeaderLen}); werr == nil {
		t.Fatal("WriteCheckpoint without RunID should fail")
	}
	got2, ok, err := LoadCheckpoint(path, "run-1")
	if err != nil || !ok || got2.StoreSeq != newCP.StoreSeq {
		t.Fatalf("checkpoint altered on validation error: (%+v, %v, %v)", got2, ok, err)
	}
}

func TestEndOffset(t *testing.T) {
	const total = 15
	const stride = 4

	path := filepath.Join(t.TempDir(), "run.swanlab")
	w := mustCreate(t, path)
	a, err := NewAsyncWriter(w, AsyncOptions{
		Queue:         32,
		IndexStride:   stride,
		FlushBytes:    1 << 20,
		FlushInterval: time.Hour,
	})
	if err != nil {
		_ = w.Close()
		t.Fatalf("NewAsyncWriter: %v", err)
	}

	// 同步写入一份独立期望偏移：第 i 条记录结束后的位置（含 1 字节边界 padding）。
	wantEnd := make([]int64, total+1)
	wantEnd[0] = HeaderLen
	for i := 1; i <= total; i++ {
		p := payloadN(10 + i*19)
		if i == 5 {
			// 第 5 条记录特意把块尾撑到只剩 4 字节，第 6 条之前必须补零。
			p = payloadN(BlockLen - int(w.Offset()%BlockLen) - HeaderLen - 4)
		}
		a.In() <- Request{StoreSeq: int64(i), Payload: p}
		if ferr := a.Flush(context.Background()); ferr != nil {
			t.Fatalf("Flush(%d): %v", i, ferr)
		}
		seq, off := a.Flushed()
		if seq != int64(i) {
			t.Fatalf("Flushed seq = %d, want %d", seq, i)
		}
		wantEnd[i] = off
	}
	defer func() { _ = a.Close() }()

	_, flushedOff := a.Flushed()
	idx := a.Index()

	t.Run("with index", func(t *testing.T) {
		for seq := int64(0); seq <= total; seq++ {
			got, oerr := EndOffset(path, idx, seq, flushedOff)
			if oerr != nil {
				t.Fatalf("EndOffset(seq=%d): %v", seq, oerr)
			}
			if got != wantEnd[seq] {
				t.Fatalf("EndOffset(seq=%d) = %d, want %d", seq, got, wantEnd[seq])
			}
		}
	})

	t.Run("without index scans from header", func(t *testing.T) {
		for _, seq := range []int64{0, 1, 5, 6, 8, total} {
			got, oerr := EndOffset(path, nil, seq, -1)
			if oerr != nil {
				t.Fatalf("EndOffset(nil, seq=%d): %v", seq, oerr)
			}
			if got != wantEnd[seq] {
				t.Fatalf("EndOffset(nil, seq=%d) = %d, want %d", seq, got, wantEnd[seq])
			}
		}
	})

	t.Run("seq beyond limit returns ErrTruncated", func(t *testing.T) {
		if _, oerr := EndOffset(path, idx, 7, wantEnd[6]); !errors.Is(oerr, ErrTruncated) {
			t.Fatalf("EndOffset beyond limit err = %v, want ErrTruncated", oerr)
		}
		if _, oerr := EndOffset(path, nil, total+1, -1); !errors.Is(oerr, ErrTruncated) {
			t.Fatalf("EndOffset beyond EOF err = %v, want ErrTruncated", oerr)
		}
	})
}

func TestEndOffsetUsesRecentCacheAndFallsBackOnMiss(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run.swanlab")
	w := mustCreate(t, path)
	a, err := NewAsyncWriter(w, AsyncOptions{
		Queue:         128,
		IndexStride:   16,
		FlushBytes:    1 << 20,
		FlushInterval: time.Hour,
	})
	if err != nil {
		_ = w.Close()
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	defer func() { _ = a.Close() }()

	// 写入 5 条记录并刷盘，验证缓存命中时不重开磁盘文件。
	for i := 1; i <= 5; i++ {
		a.In() <- Request{StoreSeq: int64(i), Payload: payloadN(20 + i)}
	}
	if ferr := a.Flush(context.Background()); ferr != nil {
		t.Fatalf("Flush: %v", ferr)
	}
	_, flushedOff := a.Flushed()
	idx := a.Index()

	cachedEnd, ok := idx.LookupEnd(3, flushedOff)
	if !ok {
		t.Fatal("expected LookupEnd(3) cache hit")
	}
	// 传不存在的路径验证 EndOffset 命中环形缓存时直接返回、不打开文件。
	missingPath := filepath.Join(t.TempDir(), "does-not-exist.swanlab")
	got, err := EndOffset(missingPath, idx, 3, flushedOff)
	if err != nil {
		t.Fatalf("EndOffset cache hit: %v", err)
	}
	if got != cachedEnd {
		t.Fatalf("EndOffset cache hit = %d, want %d", got, cachedEnd)
	}

	// limit 低于缓存结束偏移时不得误命中，应回退扫盘并按 limit 截断报错。
	if _, oerr := EndOffset(path, idx, 3, cachedEnd-1); !errors.Is(oerr, ErrTruncated) {
		t.Fatalf("EndOffset with tighter limit err = %v, want ErrTruncated", oerr)
	}

	// 写入超过 recentEndCap 条记录，验证旧序号被环表淘汰后仍能回退扫盘准确定位。
	const extra = recentEndCap + 8
	for i := 6; i <= extra; i++ {
		a.In() <- Request{StoreSeq: int64(i), Payload: payloadN(16)}
	}
	if ferr := a.Flush(context.Background()); ferr != nil {
		t.Fatalf("Flush extra: %v", ferr)
	}
	_, newFlushedOff := a.Flushed()
	if _, stillHit := idx.LookupEnd(3, newFlushedOff); stillHit {
		t.Fatal("expected seq 3 to be evicted from recent ring cache")
	}
	fallbackEnd, err := EndOffset(path, idx, 3, newFlushedOff)
	if err != nil {
		t.Fatalf("EndOffset fallback on cache miss: %v", err)
	}
	if fallbackEnd != cachedEnd {
		t.Fatalf("EndOffset fallback = %d, want %d", fallbackEnd, cachedEnd)
	}
}
