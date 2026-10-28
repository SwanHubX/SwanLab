package store

import (
	"bufio"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"
	"os"
	"path/filepath"
)

// Reader 顺序扫描日志文件，逐条返回完整记录的 payload。
//
// Scan 的错误分类是调用方判断「能否继续读」的唯一依据：
//   - io.EOF：正常读尽，或到达 OpenOffset 设置的 limit 上界；
//   - ErrTruncated：文件在记录中间结束，此前返回的记录仍然有效；
//   - ErrCorrupt：帧结构非法或 CRC 不一致，其后的数据不可信。
type Reader struct {
	fp    *os.File
	br    *bufio.Reader
	index int64
}

// Open 打开日志文件并从第一条记录开始读取，不设上界。
func Open(path string) (*Reader, error) { return OpenOffset(path, HeaderLen, -1) }

// OpenOffset 从指定字节偏移开始读取，并可选择以 limit 作为读取上界。
//
// offset 必须是某条记录的起始偏移（即 Writer.Offset 在写入前的返回值），
// 且不小于 HeaderLen。limit < 0 表示不设上界；否则读取绝不越过 limit 字节——
// 传入 AsyncWriter 发布的 flushedOffset，即可保证并发读取时看不到尚未刷盘的半条记录。
func OpenOffset(path string, offset, limit int64) (*Reader, error) {
	if offset < HeaderLen {
		return nil, fmt.Errorf("store: offset %d is before the first record (%d)", offset, HeaderLen)
	}
	fp, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	r := &Reader{fp: fp}
	if err := r.readHeader(); err != nil {
		fp.Close()
		return nil, fmt.Errorf("%s: %w", filepath.Base(path), err)
	}
	if _, err := fp.Seek(offset, io.SeekStart); err != nil {
		fp.Close()
		return nil, err
	}
	r.index = offset
	var src io.Reader = fp
	if limit >= 0 {
		if limit < offset {
			limit = offset
		}
		src = io.LimitReader(fp, limit-offset)
	}
	r.br = bufio.NewReaderSize(src, 64*1024)
	return r, nil
}

// Position 返回下一条记录起始处的文件字节偏移，用于增量扫描的进度展示。
func (r *Reader) Position() int64 { return r.index }

// Close 关闭底层文件。可重复调用。
func (r *Reader) Close() error {
	if r.fp == nil {
		return nil
	}
	err := r.fp.Close()
	r.fp = nil
	return err
}

// Scan 返回下一条完整记录的 payload。
func (r *Reader) Scan() ([]byte, error) {
	if err := r.skipPadding(); err != nil {
		return nil, err
	}
	start := r.index
	typ, payload, err := r.readChunk()
	if err != nil {
		return nil, err
	}
	switch typ {
	case chunkFull:
		return payload, nil
	case chunkFirst:
	default:
		return nil, fmt.Errorf("%w: expected FULL or FIRST chunk, got type %d at offset %d", ErrCorrupt, typ, start)
	}
	// 跨块记录：拼接 FIRST + MIDDLE* + LAST
	for {
		typ, chunk, err := r.readChunk()
		if err != nil {
			if errors.Is(err, io.EOF) {
				return nil, fmt.Errorf("%w: FIRST chunk without LAST at offset %d", ErrTruncated, start)
			}
			return nil, err
		}
		payload = append(payload, chunk...)
		switch typ {
		case chunkLast:
			return payload, nil
		case chunkMiddle:
		default:
			return nil, fmt.Errorf("%w: expected MIDDLE or LAST chunk, got type %d at offset %d", ErrCorrupt, typ, r.index)
		}
	}
}

func (r *Reader) readHeader() error {
	var h [HeaderLen]byte
	if _, err := io.ReadFull(r.fp, h[:]); err != nil {
		return fmt.Errorf("%w: %v", ErrBadHeader, err)
	}
	if ident := string(h[0:4]); ident != headerIdent {
		return fmt.Errorf("%w: ident=%q", ErrBadHeader, ident)
	}
	if magic := binary.LittleEndian.Uint16(h[4:6]); magic != headerMagic {
		return fmt.Errorf("%w: magic=0x%04X", ErrBadHeader, magic)
	}
	if v := h[6]; v != HeaderVersion {
		return fmt.Errorf("%w: version=%d (expected %d), 该文件由不兼容的 SwanLab SDK 版本生成",
			ErrBadHeader, v, HeaderVersion)
	}
	return nil
}

// skipPadding 跳过块尾的零填充。
//
// 填充只在写入下一条记录时才会补上，因此文件恰好停在记录末尾时读到 EOF 属于正常结束。
func (r *Reader) skipPadding() error {
	rest := BlockLen - r.index%BlockLen
	if rest >= HeaderLen {
		return nil
	}
	pad := make([]byte, rest)
	n, err := io.ReadFull(r.br, pad)
	if err != nil {
		if n == 0 {
			return io.EOF
		}
		return fmt.Errorf("%w: block padding at offset %d: %v", ErrTruncated, r.index, err)
	}
	for i, b := range pad {
		if b != 0 {
			return fmt.Errorf("%w: non-zero block padding at offset %d", ErrCorrupt, r.index+int64(i))
		}
	}
	r.index += int64(rest)
	return nil
}

// readChunk 读取一个 chunk 头与 payload 并校验 CRC。
// 只有连续 0 字节可读时才返回 io.EOF，其余情况都归入 ErrTruncated。
func (r *Reader) readChunk() (byte, []byte, error) {
	var h [HeaderLen]byte
	n, err := io.ReadFull(r.br, h[:])
	if err != nil {
		if n == 0 {
			return 0, nil, io.EOF
		}
		return 0, nil, fmt.Errorf("%w: chunk header at offset %d: %v", ErrTruncated, r.index, err)
	}
	start := r.index
	checksum := binary.LittleEndian.Uint32(h[0:4])
	length := int(binary.LittleEndian.Uint16(h[4:6]))
	typ := h[6]
	if typ < chunkFull || typ > chunkLast {
		return 0, nil, fmt.Errorf("%w: invalid chunk type %d at offset %d", ErrCorrupt, typ, start)
	}
	payload := make([]byte, length)
	if length > 0 {
		if _, err := io.ReadFull(r.br, payload); err != nil {
			return 0, nil, fmt.Errorf("%w: chunk payload at offset %d: %v", ErrTruncated, start, err)
		}
	}
	r.index += int64(HeaderLen + length)
	if got := crc32.Update(crcSeed[typ], crc32.IEEETable, payload); got != checksum {
		return 0, nil, fmt.Errorf("%w: crc mismatch at offset %d (got 0x%08X, want 0x%08X)",
			ErrCorrupt, start, got, checksum)
	}
	return typ, payload, nil
}
