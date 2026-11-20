package store

import (
	"errors"
	"io"
	"os"
	"path/filepath"
	"sync"
	"testing"
	"time"
)

// recordSink 把写入转发给真实文件，并记录每次 Write 的字节数；可注入延迟，用于模拟慢盘。
type recordSink struct {
	fp    *os.File
	delay time.Duration

	mu    sync.Mutex
	sizes []int
	total int
}

func (s *recordSink) Write(p []byte) (int, error) {
	if s.delay > 0 {
		time.Sleep(s.delay)
	}
	n, err := s.fp.Write(p)
	s.mu.Lock()
	s.sizes = append(s.sizes, n)
	s.total += n
	s.mu.Unlock()
	return n, err
}

func (s *recordSink) Close() error { return s.fp.Close() }

func (s *recordSink) snapshot() ([]int, int) {
	s.mu.Lock()
	defer s.mu.Unlock()
	return append([]int(nil), s.sizes...), s.total
}

// newSinkWriter 在已写好文件头的日志文件上挂一个可观测的底层写入器。
func newSinkWriter(t *testing.T, delay time.Duration) (*Writer, *recordSink, string) {
	t.Helper()
	path := filepath.Join(t.TempDir(), "run.swanlab")
	header, err := Create(path)
	if err != nil {
		t.Fatalf("Create: %v", err)
	}
	if err = header.Close(); err != nil {
		t.Fatalf("close header writer: %v", err)
	}
	fp, err := os.OpenFile(path, os.O_WRONLY|os.O_APPEND, 0o644)
	if err != nil {
		t.Fatalf("open append: %v", err)
	}
	sink := &recordSink{fp: fp, delay: delay}
	t.Cleanup(func() { _ = sink.Close() })
	return NewWriterForTest(sink, HeaderLen), sink, path
}

// scanBoundary 返回 limit 之内的完整记录条数与最后一条完整记录的结束偏移。
func scanBoundary(t *testing.T, path string, limit int64) (uint64, int64) {
	t.Helper()
	r, err := OpenOffset(path, HeaderLen, limit)
	if err != nil {
		t.Fatalf("OpenOffset(%s, %d, %d): %v", path, HeaderLen, limit, err)
	}
	defer func() { _ = r.Close() }()
	var (
		seq uint64
		end = int64(HeaderLen)
	)
	for {
		if _, err = r.Scan(); err != nil {
			if errors.Is(err, io.EOF) || errors.Is(err, ErrTruncated) {
				return seq, end
			}
			t.Fatalf("Scan(%s, limit %d): %v", path, limit, err)
		}
		seq++
		end = r.Position()
	}
}

// TestAsyncWatermarkStopsAtLastCompleteRecord 验证缓冲未排空时进度停在上一条完整
// 记录，且写入过程中没有小于整块的碎片写入。
func TestAsyncWatermarkStopsAtLastCompleteRecord(t *testing.T) {
	w, sink, path := newSinkWriter(t, 0)
	a, err := NewAsyncWriter(w, AsyncOptions{Queue: 64, FlushInterval: time.Hour})
	if err != nil {
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	defer func() { _ = a.Close() }()

	const n = 1200
	for i := 0; i < n; i++ {
		a.In() <- Request{StoreSeq: uint64(i + 1), Payload: payloadN(50)}
	}
	// 只写完、不刷盘：缓冲里仍有数据，进度只能停在上一条完整记录。
	if err = a.sync(false); err != nil {
		t.Fatalf("sync(false): %v", err)
	}
	written := a.w.written()
	if a.w.Offset() <= written {
		t.Fatalf("buffer should still hold bytes: offset=%d written=%d", a.w.Offset(), written)
	}

	seq, off := a.Flushed()
	wantSeq, wantOff := scanBoundary(t, path, written)
	if seq != wantSeq || off != wantOff {
		t.Fatalf("Flushed = (%d, %d), want (%d, %d)", seq, off, wantSeq, wantOff)
	}
	if seq >= n {
		t.Fatalf("flushed seq = %d, want progress to lag behind %d records", seq, n)
	}

	sizes, total := sink.snapshot()
	if got := int64(total) + HeaderLen; got != written {
		t.Fatalf("observed %d bytes, writer reports %d", got, written)
	}
	for i, size := range sizes {
		if size != BlockLen {
			t.Fatalf("write %d carried %d bytes, want a full block (%d)", i, size, BlockLen)
		}
	}
}

// TestAsyncWatermarkHidesPartialCrossBlockRecord 验证跨块大记录只写了一部分时，
// 进度不会推到半条记录上。
func TestAsyncWatermarkHidesPartialCrossBlockRecord(t *testing.T) {
	w, _, path := newSinkWriter(t, 0)
	a, err := NewAsyncWriter(w, AsyncOptions{Queue: 64, FlushInterval: time.Hour})
	if err != nil {
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	defer func() { _ = a.Close() }()

	// 先写满一个块，再追加一条跨块记录，使其末尾字节仍留在缓冲里。
	for i := 0; i < 600; i++ {
		a.In() <- Request{StoreSeq: uint64(i + 1), Payload: payloadN(50)}
	}
	a.In() <- Request{StoreSeq: 601, Payload: payloadN(40 * 1024)}
	if err = a.sync(false); err != nil {
		t.Fatalf("sync(false): %v", err)
	}

	seq, off := a.Flushed()
	wantSeq, wantOff := scanBoundary(t, path, a.w.written())
	if seq != wantSeq || off != wantOff {
		t.Fatalf("Flushed = (%d, %d), want (%d, %d)", seq, off, wantSeq, wantOff)
	}
	if seq >= 601 {
		t.Fatalf("flushed seq = %d, want the last record to stay buffered", seq)
	}
	if got := readLimit(t, path, HeaderLen, off); len(got) != int(seq) {
		t.Fatalf("read %d records below flushed offset, want %d", len(got), seq)
	}
}

// TestAsyncSlowDiskWatermarkStaysWithinWritten 验证慢盘下进度始终不超过真实写出的
// 位置、只对应完整记录，且不会回退。
func TestAsyncSlowDiskWatermarkStaysWithinWritten(t *testing.T) {
	w, _, path := newSinkWriter(t, time.Millisecond)
	a, err := NewAsyncWriter(w, AsyncOptions{Queue: 8, FlushInterval: time.Hour})
	if err != nil {
		t.Fatalf("NewAsyncWriter: %v", err)
	}

	const n = 300
	done := make(chan error, 1)
	go func() {
		for i := 0; i < n; i++ {
			a.In() <- Request{StoreSeq: uint64(i + 1), Payload: payloadN(80)}
		}
		done <- a.Close()
	}()

	var (
		lastSeq uint64
		lastOff = int64(HeaderLen)
	)
	for closed := false; !closed; {
		select {
		case err = <-done:
			if err != nil {
				t.Fatalf("Close: %v", err)
			}
			closed = true
		default:
		}
		seq, off := a.Flushed()
		if seq < lastSeq || off < lastOff {
			t.Fatalf("watermark regressed: (%d, %d) after (%d, %d)", seq, off, lastSeq, lastOff)
		}
		lastSeq, lastOff = seq, off
		fi, statErr := os.Stat(path)
		if statErr != nil {
			t.Fatalf("Stat: %v", statErr)
		}
		if off > fi.Size() {
			t.Fatalf("flushed offset %d beyond file size %d", off, fi.Size())
		}
		if got, end := scanBoundary(t, path, off); got != seq || end != off {
			t.Fatalf("watermark (%d, %d) covers (%d, %d)", seq, off, got, end)
		}
	}

	fi, err := os.Stat(path)
	if err != nil {
		t.Fatalf("Stat: %v", err)
	}
	if seq, off := a.Flushed(); seq != n || off != fi.Size() {
		t.Fatalf("final Flushed = (%d, %d), want (%d, %d)", seq, off, n, fi.Size())
	}
}

// TestAsyncTickerSkipsCyclesWithNaturalWrites 验证定时器只在本周期没有自然写出时才
// 刷盘。直接驱动写协程的状态，避开对定时器到点时机的依赖。
func TestAsyncTickerSkipsCyclesWithNaturalWrites(t *testing.T) {
	w, sink, _ := newSinkWriter(t, 0)
	a, err := NewAsyncWriter(w, AsyncOptions{Queue: 64, FlushInterval: time.Hour})
	if err != nil {
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	defer func() { _ = a.Close() }()

	feed := func(from, to int) {
		for i := from; i < to; i++ {
			a.In() <- Request{StoreSeq: uint64(i + 1), Payload: payloadN(200)}
		}
		if syncErr := a.sync(false); syncErr != nil {
			t.Fatalf("sync(false): %v", syncErr)
		}
	}

	// 先写到发生一次自然写出，并假装上一个定时周期已经看到它。
	feed(0, 400)
	written := a.w.written()
	if written == HeaderLen {
		t.Fatalf("no natural write happened: written=%d", written)
	}
	if a.w.Offset() <= written {
		t.Fatalf("buffer should still hold bytes: offset=%d written=%d", a.w.Offset(), written)
	}
	st := writeState{seq: 400, lastWritten: written, lastFlush: written}

	// 本周期内又发生一次自然写出：这次定时检查不得插队刷盘。
	feed(400, 800)
	advancer := a.w.written()
	if advancer <= written {
		t.Fatalf("written = %d, want progress beyond %d", advancer, written)
	}
	if a.w.Offset() <= advancer {
		t.Fatalf("buffer should still hold bytes: offset=%d written=%d", a.w.Offset(), advancer)
	}
	before := len(sinkSizes(t, sink))
	a.tick(&st)
	if got := len(sinkSizes(t, sink)); got != before {
		t.Fatalf("ticker flushed %d extra writes after a natural write in the same cycle", got-before)
	}
	if st.lastWritten != advancer {
		t.Fatalf("lastWritten = %d, want %d", st.lastWritten, advancer)
	}

	// 本周期没有任何自然写出、但缓冲非空：必须刷盘。
	// st.seq 与写协程一致地推进到最后一条记录，刷盘后进度才能覆盖它。
	st.seq = 800
	a.tick(&st)
	if got := len(sinkSizes(t, sink)); got != before+1 {
		t.Fatalf("ticker wrote %d times, want exactly 1 fallback flush", got-before)
	}
	if st.lastWritten != a.w.Offset() {
		t.Fatalf("lastWritten = %d, want %d after fallback flush", st.lastWritten, a.w.Offset())
	}
	if seq, off := a.Flushed(); seq != 800 || off != a.w.Offset() {
		t.Fatalf("Flushed = (%d, %d), want (800, %d)", seq, off, a.w.Offset())
	}
}

func sinkSizes(t *testing.T, sink *recordSink) []int {
	t.Helper()
	sizes, _ := sink.snapshot()
	return sizes
}

// TestAsyncTickerFlushesIdleBuffer 验证低吞吐停写后定时器会把缓冲写入。
func TestAsyncTickerFlushesIdleBuffer(t *testing.T) {
	w, _, path := newSinkWriter(t, 0)
	a, err := NewAsyncWriter(w, AsyncOptions{Queue: 8, FlushInterval: 5 * time.Millisecond})
	if err != nil {
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	defer func() { _ = a.Close() }()

	a.In() <- Request{StoreSeq: 1, Payload: payloadN(20)}
	deadline := time.Now().Add(2 * time.Second)
	for {
		seq, off := a.Flushed()
		if seq == 1 {
			fi, statErr := os.Stat(path)
			if statErr != nil {
				t.Fatalf("Stat: %v", statErr)
			}
			if off != fi.Size() {
				t.Fatalf("flushed offset = %d, file size = %d", off, fi.Size())
			}
			return
		}
		if time.Now().After(deadline) {
			t.Fatal("timed flush did not advance progress")
		}
		time.Sleep(time.Millisecond)
	}
}

// shortSink 每次最多只接受 n 字节，用于模拟底层短写。
type shortSink struct {
	n int
}

func (s *shortSink) Write(p []byte) (int, error) {
	if s.n >= len(p) {
		return len(p), nil
	}
	return s.n, io.ErrShortWrite
}

func (s *shortSink) Close() error { return nil }

func TestCountingWriterCountsActualBytes(t *testing.T) {
	cw := &countingWriter{fp: &shortSink{n: 3}}
	n, err := cw.Write([]byte("0123456789"))
	if n != 3 || !errors.Is(err, io.ErrShortWrite) {
		t.Fatalf("Write = (%d, %v), want (3, ErrShortWrite)", n, err)
	}
	if got := cw.written.Load(); got != 3 {
		t.Fatalf("written = %d, want 3", got)
	}
}

func TestAsyncShortWriteDoesNotAdvanceWatermark(t *testing.T) {
	sink := &shortSink{n: 10}
	a, err := NewAsyncWriter(NewWriterForTest(sink, HeaderLen), AsyncOptions{Queue: 4, FlushInterval: time.Hour})
	if err != nil {
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	defer func() { _ = a.Close() }()

	a.In() <- Request{StoreSeq: 1, Payload: payloadN(2000)}
	if err = a.FlushForTest(); !errors.Is(err, io.ErrShortWrite) {
		t.Fatalf("FlushForTest err = %v, want ErrShortWrite", err)
	}
	if seq, off := a.Flushed(); seq != 0 || off != HeaderLen {
		t.Fatalf("Flushed = (%d, %d) after short write, want (0, %d)", seq, off, HeaderLen)
	}
	if got := a.w.written(); got != HeaderLen+10 {
		t.Fatalf("written = %d, want %d", got, HeaderLen+10)
	}
}
