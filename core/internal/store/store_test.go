package store

import (
	"bytes"
	"encoding/binary"
	"errors"
	"hash/crc32"
	"io"
	"os"
	"path/filepath"
	"strings"
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

// writeFile 写出一组记录并返回各条记录的起始偏移（含文件头偏移）。
func writeFile(t *testing.T, payloads ...[]byte) (string, []int64) {
	t.Helper()
	path := filepath.Join(t.TempDir(), "run.swanlab")
	w := mustCreate(t, path)
	offsets := make([]int64, len(payloads))
	for i, p := range payloads {
		offsets[i] = w.Offset()
		if err := w.Write(p); err != nil {
			t.Fatalf("Write record %d: %v", i, err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return path, offsets
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
		want := make([][]byte, 0, len(sizes))
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

	want := make([]byte, 0, HeaderLen+HeaderLen+len(payload))
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
	wantEnd := firstEnd + (BlockLen - firstEnd%BlockLen) + HeaderLen + int64(len(second))
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

func TestOpenOffsetAndLimit(t *testing.T) {
	payloads := [][]byte{[]byte("aaa"), []byte("bbbb"), []byte("ccccc")}
	path, off := writeFile(t, payloads...)

	// 从第二条记录起读，上界正好是第二条记录的结束偏移
	got := readLimit(t, path, off[1], off[2])
	if len(got) != 1 || !bytes.Equal(got[0], payloads[1]) {
		t.Fatalf("got %d records, want only record 1", len(got))
	}

	// limit 小于 offset 时视为空窗口
	if empty := readLimit(t, path, off[2], off[1]); len(empty) != 0 {
		t.Fatalf("empty window returned %d records", len(empty))
	}

	// 无上界时可读到文件末尾
	got = readLimit(t, path, off[1], -1)
	if len(got) != 2 || !bytes.Equal(got[1], payloads[2]) {
		t.Fatalf("open-ended read from offset %d got %d records", off[1], len(got))
	}
}

func TestLimitInsideRecord(t *testing.T) {
	payloads := [][]byte{[]byte("first-record"), []byte("second-record")}
	path, off := writeFile(t, payloads...)

	// 上界落在第二条记录的 payload 中间：不得返回半条记录
	r, err := OpenOffset(path, off[1], off[1]+HeaderLen+4)
	if err != nil {
		t.Fatalf("OpenOffset: %v", err)
	}
	defer r.Close()
	if p, err := r.Scan(); !errors.Is(err, ErrTruncated) {
		t.Fatalf("Scan = %q, %v; want ErrTruncated", p, err)
	}
}

func TestOpenOffsetBeforeFirstRecord(t *testing.T) {
	path, _ := writeFile(t, []byte("x"))
	if _, err := OpenOffset(path, HeaderLen-1, -1); err == nil {
		t.Fatal("OpenOffset with offset < HeaderLen should fail")
	}
}

func TestOpenMissingFile(t *testing.T) {
	if _, err := Open(filepath.Join(t.TempDir(), "nope.swanlab")); err == nil {
		t.Fatal("Open of missing file should fail")
	}
}

func TestBadHeader(t *testing.T) {
	cases := []struct {
		name   string
		header []byte
		want   string
	}{
		{"empty", nil, ""},
		{"ident", []byte{'X', 'S', 'W', 'L', 0xD6, 0xE1, HeaderVersion}, "ident"},
		{"magic", []byte{':', 'S', 'W', 'L', 0x00, 0xE1, HeaderVersion}, "magic"},
		{"version", []byte{':', 'S', 'W', 'L', 0xD6, 0xE1, 1}, "invalid run version: 1"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "run.swanlab")
			if err := os.WriteFile(path, c.header, 0o644); err != nil {
				t.Fatalf("WriteFile: %v", err)
			}
			_, err := Open(path)
			if !errors.Is(err, ErrBadHeader) {
				t.Fatalf("Open error = %v, want ErrBadHeader", err)
			}
			if c.want != "" && !strings.Contains(err.Error(), c.want) {
				t.Fatalf("error %q does not mention %q", err, c.want)
			}
		})
	}
}

func TestReaderTruncated(t *testing.T) {
	payloads := [][]byte{[]byte("first-record"), []byte("second-record"), bytes.Repeat([]byte("z"), 100)}
	path, off := writeFile(t, payloads...)
	fi, err := os.Stat(path)
	if err != nil {
		t.Fatalf("Stat: %v", err)
	}

	cases := []struct {
		name     string
		size     int64
		complete int
	}{
		{"mid payload of last", fi.Size() - 3, 2},
		{"mid chunk header of last", off[2] + 3, 2},
		{"mid payload of first", off[0] + HeaderLen + 4, 0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			truncated := filepath.Join(t.TempDir(), "trunc.swanlab")
			data, err := os.ReadFile(path)
			if err != nil {
				t.Fatalf("ReadFile: %v", err)
			}
			if err = os.WriteFile(truncated, data[:c.size], 0o644); err != nil {
				t.Fatalf("WriteFile: %v", err)
			}
			r, err := Open(truncated)
			if err != nil {
				t.Fatalf("Open: %v", err)
			}
			defer r.Close()

			complete := 0
			var lastErr error
			for {
				if _, err := r.Scan(); err != nil {
					lastErr = err
					break
				}
				complete++
			}
			if !errors.Is(lastErr, ErrTruncated) {
				t.Fatalf("scan stop error = %v, want ErrTruncated", lastErr)
			}
			if complete != c.complete {
				t.Fatalf("stopped after %d complete records, want %d", complete, c.complete)
			}
		})
	}

	// 恰好截在记录边界：属于正常读完
	exact := filepath.Join(t.TempDir(), "exact.swanlab")
	data, _ := os.ReadFile(path)
	if err := os.WriteFile(exact, data[:off[2]], 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	if got := readAll(t, exact); len(got) != 2 {
		t.Fatalf("boundary truncation returned %d records, want 2", len(got))
	}
}

func TestReaderCorrupt(t *testing.T) {
	payloads := [][]byte{[]byte("first-record"), []byte("second-record")}
	path, off := writeFile(t, payloads...)
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}

	corrupt := func(name string, mutate func([]byte)) {
		t.Run(name, func(t *testing.T) {
			b := append([]byte(nil), data...)
			mutate(b)
			p := filepath.Join(t.TempDir(), "corrupt.swanlab")
			if err := os.WriteFile(p, b, 0o644); err != nil {
				t.Fatalf("WriteFile: %v", err)
			}
			r, err := Open(p)
			if err != nil {
				t.Fatalf("Open: %v", err)
			}
			defer r.Close()
			if _, err := r.Scan(); err != nil {
				t.Fatalf("first record should still be readable: %v", err)
			}
			if _, err := r.Scan(); !errors.Is(err, ErrCorrupt) {
				t.Fatalf("second Scan error = %v, want ErrCorrupt", err)
			}
		})
	}

	// payload 字节翻转 → CRC 不一致
	corrupt("checksum", func(b []byte) { b[off[1]+HeaderLen] ^= 0xFF })
	// chunk 类型非法
	corrupt("chunk type", func(b []byte) { b[off[1]+HeaderLen-1] = 9 })
}

func TestReaderNonZeroPadding(t *testing.T) {
	path, _ := writeFile(t, payloadN(BlockLen-HeaderLen*2-6), payloadN(10))
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	// 第一条记录结束于 32762，其后 6 字节为零填充
	if data[BlockLen-1] != 0 {
		t.Fatalf("expected zero padding at block tail, got %#x", data[BlockLen-1])
	}
	data[BlockLen-1] = 0x01
	p := filepath.Join(t.TempDir(), "pad.swanlab")
	if err = os.WriteFile(p, data, 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	r, err := Open(p)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer r.Close()
	if _, err := r.Scan(); err != nil {
		t.Fatalf("first record: %v", err)
	}
	if _, err := r.Scan(); !errors.Is(err, ErrCorrupt) {
		t.Fatalf("Scan after bad padding = %v, want ErrCorrupt", err)
	}
}

func TestReaderPosition(t *testing.T) {
	payloads := [][]byte{[]byte("aaa"), []byte("bbbb")}
	path, off := writeFile(t, payloads...)
	r, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer r.Close()
	if got := r.Position(); got != HeaderLen {
		t.Fatalf("initial position = %d, want %d", got, HeaderLen)
	}
	if _, err := r.Scan(); err != nil {
		t.Fatalf("Scan: %v", err)
	}
	if got := r.Position(); got != off[1] {
		t.Fatalf("position after first = %d, want %d", got, off[1])
	}
}
