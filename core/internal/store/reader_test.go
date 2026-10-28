package store

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

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
