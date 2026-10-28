package store

import (
	"bytes"
	"errors"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"testing"
	"time"
)

func newTestAsync(t *testing.T, opts AsyncOptions) (*AsyncWriter, string) {
	t.Helper()
	path := filepath.Join(t.TempDir(), "run.swanlab")
	w := mustCreate(t, path)
	a, err := NewAsyncWriter(w, opts)
	if err != nil {
		w.Close()
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	return a, path
}

func anyPayload(i int) []byte {
	// 尺寸在 1..~300 字节之间变化，覆盖单块内的各种切分
	return payloadN(1 + i*37%300)
}

func TestAsyncWriteAndClose(t *testing.T) {
	a, path := newTestAsync(t, AsyncOptions{Queue: 64})

	const n = 200
	for i := 0; i < n; i++ {
		a.In() <- Request{StoreSeq: uint64(i + 1), Payload: anyPayload(i)}
	}
	if err := a.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	got := readAll(t, path)
	if len(got) != n {
		t.Fatalf("read %d records, want %d", len(got), n)
	}
	for i := range got {
		if !bytes.Equal(got[i], anyPayload(i)) {
			t.Fatalf("record %d mismatch", i)
		}
	}

	fi, err := os.Stat(path)
	if err != nil {
		t.Fatalf("Stat: %v", err)
	}
	seq, off := a.Flushed()
	if seq != n {
		t.Fatalf("flushed seq = %d, want %d", seq, n)
	}
	if off != fi.Size() {
		t.Fatalf("flushed offset = %d, file size = %d", off, fi.Size())
	}

	// 幂等：再次 Close 只等待同一结果
	if err := a.Close(); err != nil {
		t.Fatalf("second Close: %v", err)
	}
	if err := a.FlushForTest(); err == nil {
		t.Fatal("FlushForTest after Close should fail")
	}
}

func TestAsyncPublishesOnlyAfterFlush(t *testing.T) {
	a, path := newTestAsync(t, AsyncOptions{Queue: 8, FlushInterval: time.Hour})
	defer a.Close()

	if seq, off := a.Flushed(); seq != 0 || off != HeaderLen {
		t.Fatalf("initial flushed = (%d, %d), want (0, %d)", seq, off, HeaderLen)
	}

	for i := 0; i < 3; i++ {
		a.In() <- Request{StoreSeq: uint64(i + 1), Payload: anyPayload(i)}
	}
	// 只写完、不刷盘：记录还在用户态缓冲里，进度不得前移
	if err := a.sync(false); err != nil {
		t.Fatalf("sync(false): %v", err)
	}
	if seq, _ := a.Flushed(); seq != 0 {
		t.Fatalf("flushed seq = %d before flush, want 0", seq)
	}
	if fi, err := os.Stat(path); err != nil || fi.Size() != HeaderLen {
		t.Fatalf("file size before flush = %v, want %d (err %v)", fi.Size(), HeaderLen, err)
	}

	if err := a.FlushForTest(); err != nil {
		t.Fatalf("FlushForTest: %v", err)
	}
	fi, err := os.Stat(path)
	if err != nil {
		t.Fatalf("Stat: %v", err)
	}
	if seq, off := a.Flushed(); seq != 3 || off != fi.Size() {
		t.Fatalf("flushed = (%d, %d), want (3, %d)", seq, off, fi.Size())
	}
}

func TestAsyncQuantitativeFlush(t *testing.T) {
	// 阈值 200 字节：第一条不够，第二条越过后立即刷盘
	a, _ := newTestAsync(t, AsyncOptions{Queue: 8, FlushBytes: 200, FlushInterval: time.Hour})
	defer a.Close()

	a.In() <- Request{StoreSeq: 1, Payload: payloadN(50)}
	if err := a.sync(false); err != nil {
		t.Fatalf("sync(false): %v", err)
	}
	if seq, _ := a.Flushed(); seq != 0 {
		t.Fatalf("flushed seq = %d after 57 bytes, want 0", seq)
	}

	a.In() <- Request{StoreSeq: 2, Payload: payloadN(200)}
	if err := a.sync(false); err != nil {
		t.Fatalf("sync(false): %v", err)
	}
	if seq, _ := a.Flushed(); seq != 2 {
		t.Fatalf("flushed seq = %d after crossing threshold, want 2", seq)
	}
}

func TestAsyncTimedFlush(t *testing.T) {
	a, _ := newTestAsync(t, AsyncOptions{Queue: 8, FlushInterval: 10 * time.Millisecond})
	defer a.Close()

	a.In() <- Request{StoreSeq: 1, Payload: payloadN(20)}
	deadline := time.Now().Add(2 * time.Second)
	for {
		if seq, _ := a.Flushed(); seq == 1 {
			return
		}
		if time.Now().After(deadline) {
			t.Fatal("timed flush did not advance progress")
		}
		time.Sleep(time.Millisecond)
	}
}

func TestAsyncErrorCallback(t *testing.T) {
	var calls int
	a, _ := newTestAsync(t, AsyncOptions{Queue: 8, OnError: func(error) { calls++ }})

	// 提前关闭底层文件，模拟运行期 I/O 失败
	if err := a.w.Close(); err != nil {
		t.Fatalf("close underlying writer: %v", err)
	}
	a.In() <- Request{StoreSeq: 1, Payload: payloadN(10)}
	if err := a.FlushForTest(); err == nil {
		t.Fatal("FlushForTest after write failure should report error")
	}
	first := a.Err()
	if first == nil {
		t.Fatal("Err() = nil after write failure")
	}
	if calls != 1 {
		t.Fatalf("OnError called %d times, want 1", calls)
	}

	// 降级后继续丢弃写入，不阻断上游；错误保持首次值
	a.In() <- Request{StoreSeq: 2, Payload: payloadN(10)}
	if err := a.sync(false); err != nil {
		t.Fatalf("sync(false) after degrade: %v", err)
	}
	if a.Err() != first {
		t.Fatalf("Err() changed to %v", a.Err())
	}
	if calls != 1 {
		t.Fatalf("OnError called %d times after degrade, want 1", calls)
	}
	if seq, off := a.Flushed(); seq != 0 || off != HeaderLen {
		t.Fatalf("flushed = (%d, %d) after degrade, want (0, %d)", seq, off, HeaderLen)
	}
	if err := a.Close(); err != first {
		t.Fatalf("Close error = %v, want %v", err, first)
	}
}

func TestAsyncConcurrentReaderBoundedByFlushed(t *testing.T) {
	a, path := newTestAsync(t, AsyncOptions{Queue: 16, FlushInterval: 2 * time.Millisecond})

	const n = 1200
	writerDone := make(chan error, 1)
	go func() {
		for i := 0; i < n; i++ {
			a.In() <- Request{StoreSeq: uint64(i + 1), Payload: anyPayload(i)}
		}
		writerDone <- a.Close()
	}()

	for {
		seq, off := a.Flushed()
		if seq > 0 {
			r, err := OpenOffset(path, HeaderLen, off)
			if err != nil {
				t.Fatalf("OpenOffset: %v", err)
			}
			var count uint64
			for {
				p, err := r.Scan()
				if errors.Is(err, io.EOF) {
					break
				}
				if err != nil {
					r.Close()
					t.Fatalf("Scan while writing (flushed seq %d, offset %d): %v", seq, off, err)
				}
				if !bytes.Equal(p, anyPayload(int(count))) {
					r.Close()
					t.Fatalf("record %d mismatch while writing", count)
				}
				count++
			}
			r.Close()
			if count != seq {
				t.Fatalf("reader saw %d records, flushed seq %d", count, seq)
			}
		}

		select {
		case err := <-writerDone:
			if err != nil {
				t.Fatalf("Close: %v", err)
			}
			fi, err := os.Stat(path)
			if err != nil {
				t.Fatalf("Stat: %v", err)
			}
			seq, off := a.Flushed()
			if seq != n || off != fi.Size() {
				t.Fatalf("final flushed = (%d, %d), want (%d, %d)", seq, off, n, fi.Size())
			}
			if got := readAll(t, path); len(got) != n {
				t.Fatalf("final read %d records, want %d", len(got), n)
			}
			return
		default:
			runtime.Gosched()
		}
	}
}
