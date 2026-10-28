package store

import (
	"bytes"
	"encoding/binary"
	"errors"
	"hash/crc32"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// ---------- helpers ----------

func mustCreate(t *testing.T, path string) *Writer {
	t.Helper()
	w, err := Create(path)
	if err != nil {
		t.Fatalf("Create(%s): %v", path, err)
	}
	return w
}

func payloadN(n int) []byte {
	p := make([]byte, n)
	for i := range p {
		p[i] = byte(i*7 + n)
	}
	return p
}

// readLimit 从 offset 读到 limit，正常读尽返回 nil error。
func readLimit(t *testing.T, path string, offset, limit int64) [][]byte {
	t.Helper()
	r, err := OpenOffset(path, offset, limit)
	if err != nil {
		t.Fatalf("OpenOffset(%s, %d, %d): %v", path, offset, limit, err)
	}
	defer r.Close()
	var out [][]byte
	for {
		p, err := r.Scan()
		if errors.Is(err, io.EOF) {
			return out
		}
		if err != nil {
			t.Fatalf("Scan: %v", err)
		}
		out = append(out, p)
	}
}

func readAll(t *testing.T, path string) [][]byte {
	t.Helper()
	return readLimit(t, path, HeaderLen, -1)
}

// ---------- 读写往返 ----------

func TestRoundTrip(t *testing.T) {
	// 覆盖：空记录、单块、块剩余空间恰好等于/小于 chunk 头、跨两块、跨三块
	sizes := []int{0, 1, 17, 32740, 32747, 32748, 32760, 32761, 32762, 65522, 65523, 100000, 200000}

	t.Run("each file", func(t *testing.T) {
		// 每个尺寸独占一个文件：首条记录就要处理块尾/跨块
		for _, n := range sizes {
			path := filepath.Join(t.TempDir(), "run.swanlab")
			w := mustCreate(t, path)
			if err := w.Write(payloadN(n)); err != nil {
				t.Fatalf("size %d: Write: %v", n, err)
			}
			if err := w.Close(); err != nil {
				t.Fatalf("size %d: Close: %v", n, err)
			}
			got := readAll(t, path)
			if len(got) != 1 || !bytes.Equal(got[0], payloadN(n)) {
				t.Fatalf("size %d: round trip mismatch (records=%d)", n, len(got))
			}
			if fi, err := os.Stat(path); err != nil || fi.Size() != w.Offset() {
				t.Fatalf("size %d: file size %v, offset %d, err %v", n, fi, w.Offset(), err)
			}
		}
	})

	t.Run("one file", func(t *testing.T) {
		// 同一文件连续写：覆盖记录之间的块尾填充与跨块拼接
		path := filepath.Join(t.TempDir(), "run.swanlab")
		w := mustCreate(t, path)
		var want [][]byte
		for _, n := range sizes {
			p := payloadN(n)
			want = append(want, p)
			if err := w.Write(p); err != nil {
				t.Fatalf("size %d: Write: %v", n, err)
			}
		}
		if err := w.Close(); err != nil {
			t.Fatalf("Close: %v", err)
		}
		got := readAll(t, path)
		if len(got) != len(want) {
			t.Fatalf("got %d records, want %d", len(got), len(want))
		}
		for i := range want {
			if !bytes.Equal(got[i], want[i]) {
				t.Fatalf("record %d mismatch (len %d vs %d)", i, len(got[i]), len(want[i]))
			}
		}
	})
}

// ---------- 字节布局 ----------

func TestFileLayout(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run.swanlab")
	payload := []byte("hello")
	w := mustCreate(t, path)
	if err := w.Write(payload); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	got, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}

	var want []byte
	want = append(want, ':', 'S', 'W', 'L', 0xD6, 0xE1, HeaderVersion)
	var hdr [HeaderLen]byte
	// CRC 覆盖「类型字节 + payload」，等价于以 crc32(类型字节) 为预置基值
	binary.LittleEndian.PutUint32(hdr[0:4], crc32.ChecksumIEEE(append([]byte{chunkFull}, payload...)))
	binary.LittleEndian.PutUint16(hdr[4:6], uint16(len(payload)))
	hdr[6] = chunkFull
	want = append(want, hdr[:]...)
	want = append(want, payload...)

	if !bytes.Equal(got, want) {
		t.Fatalf("layout mismatch\n got: %x\nwant: %x", got, want)
	}
}

// ---------- 偏移记账 ----------

func TestWriterOffset(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run.swanlab")
	w := mustCreate(t, path)
	if off := w.Offset(); off != HeaderLen {
		t.Fatalf("initial offset = %d, want %d", off, HeaderLen)
	}

	// 首条记录把块尾撑到只剩 6 字节，下一条记录前需要补 6 字节零
	first := payloadN(BlockLen - HeaderLen*2 - 6)
	if err := w.Write(first); err != nil {
		t.Fatalf("Write first: %v", err)
	}
	firstEnd := int64(HeaderLen*2 + len(first))
	if off := w.Offset(); off != firstEnd {
		t.Fatalf("offset after first = %d, want %d", off, firstEnd)
	}

	second := payloadN(10)
	if err := w.Write(second); err != nil {
		t.Fatalf("Write second: %v", err)
	}
	wantEnd := firstEnd + int64(BlockLen-firstEnd%BlockLen) + HeaderLen + int64(len(second))
	if off := w.Offset(); off != wantEnd {
		t.Fatalf("offset after second = %d, want %d (padding not accounted)", off, wantEnd)
	}

	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if fi, err := os.Stat(path); err != nil || fi.Size() != wantEnd {
		t.Fatalf("file size %v, want %d (err %v)", fi, wantEnd, err)
	}
	got := readAll(t, path)
	if len(got) != 2 || !bytes.Equal(got[0], first) || !bytes.Equal(got[1], second) {
		t.Fatalf("round trip after padding mismatch (records=%d)", len(got))
	}
}

// ---------- Writer 生命周期 ----------

func TestCreateExclusive(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run.swanlab")
	w := mustCreate(t, path)
	defer w.Close()

	if _, err := Create(path); !errors.Is(err, os.ErrExist) {
		t.Fatalf("second Create error = %v, want ErrExist", err)
	}
}

func TestWriterClosed(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run.swanlab")
	w := mustCreate(t, path)
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("second Close: %v", err)
	}
	if err := w.Write([]byte("x")); !errors.Is(err, ErrClosed) {
		t.Fatalf("Write after Close = %v, want ErrClosed", err)
	}
	if err := w.Flush(); !errors.Is(err, ErrClosed) {
		t.Fatalf("Flush after Close = %v, want ErrClosed", err)
	}
}
