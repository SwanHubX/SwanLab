// Package store 实现 SwanLab 运行日志（run-*.swanlab）的本地持久化读写。
//
// 文件格式与 Python 侧 swanlab/sdk/internal/core_python/store 完全一致：
//   - 文件头 7 字节：":SWL" + uint16 小端魔数 0xE1D6 + uint8 版本号 2
//   - 其后为 LevelDB log 分帧：整个字节流按 32 KiB 分块，块内紧凑排列 chunk
//   - chunk 头 7 字节：CRC32(4B) + payload 长度(2B) + chunk 类型(1B)
//   - chunk 类型 1=FULL、2=FIRST、3=MIDDLE、4=LAST；跨块记录拆成 FIRST + MIDDLE* + LAST
//   - 块尾剩余不足 7 字节时以 0 填充，对齐到下一块
//   - CRC 校验范围是「类型字节 + payload」，即以 crc32(类型字节) 为预置基值
//
// payload 是序列化后的 swanlab.record.v1.Record。本包只负责字节搬运，
// 不 import proto，也不解释记录语义。
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
	"sync"
	"sync/atomic"
)

const (
	// HeaderLen 是文件头与 chunk 头的长度，两者恰好都是 7 字节。
	HeaderLen = 7
	// BlockLen 是物理块大小。
	BlockLen = 32768
	// HeaderVersion 是当前文件格式版本号，与 Python 侧保持一致。
	HeaderVersion = 2

	// dataLen 是单个 chunk 能承载的最大 payload 长度。
	dataLen = BlockLen - HeaderLen

	headerIdent = ":SWL"
	headerMagic = 0xE1D6

	chunkFull   byte = 1
	chunkFirst  byte = 2
	chunkMiddle byte = 3
	chunkLast   byte = 4
)

var (
	// ErrBadHeader 表示文件头不符合 :SWL 格式，或版本号不兼容。
	ErrBadHeader = errors.New("store: invalid file header")
	// ErrCorrupt 表示帧结构非法或 CRC 校验失败，其后的数据已不可信。
	ErrCorrupt = errors.New("store: corrupt record")
	// ErrTruncated 表示文件在记录中间结束（写入方尚未写完或进程异常退出）。
	ErrTruncated = errors.New("store: truncated record")
	// ErrClosed 表示操作发生在已关闭的 Writer 或 AsyncWriter 上。
	ErrClosed = errors.New("store: closed")
)

// crcSeed[t] 是类型字节 t 的 CRC32，作为该 chunk 校验和的预置基值。
var crcSeed = func() (s [5]uint32) {
	for t := chunkFull; t <= chunkLast; t++ {
		s[t] = crc32.ChecksumIEEE([]byte{t})
	}
	return s
}()

// Writer 以排他方式创建日志文件，并按 LevelDB log 格式顺序追加记录。
//
// Writer 只负责字节落盘与偏移记账，不解释 payload 语义，也不决定何时刷盘；
// 刷盘节奏由 AsyncWriter 控制。写入方法并发安全，但块内偏移只在单写者下才有意义，
// 因此仍应只由一个 goroutine 调用 Write。任意一次写入或刷盘失败后，
// Writer 的缓冲与文件可能已不一致，应视为不可用并整体关闭。
type Writer struct {
	mu     sync.Mutex
	fp     *os.File
	buf    *bufio.Writer
	closed bool
	index  atomic.Int64 // 当前文件末尾偏移，含文件头与块尾填充
}

// Create 排他创建日志文件并写入 7 字节文件头。
//
// 文件已存在时返回的错误满足 errors.Is(err, fs.ErrExist)；调用方据此识别
// 同一 run_dir 被重复打开或 Resume 冲突。
func Create(path string) (*Writer, error) {
	fp, err := os.OpenFile(path, os.O_WRONLY|os.O_CREATE|os.O_EXCL, 0o644)
	if err != nil {
		return nil, err
	}
	w := &Writer{fp: fp, buf: bufio.NewWriterSize(fp, BlockLen)}
	w.index.Store(HeaderLen)
	if err := w.writeHeader(); err != nil {
		fp.Close()
		return nil, err
	}
	// 文件头立即刷出：AsyncWriter 发布的初始 (flushedStoreSeq, flushedOffset)
	// 以「文件头对读取方可见」为前提。
	if err := w.buf.Flush(); err != nil {
		fp.Close()
		return nil, err
	}
	return w, nil
}

func (w *Writer) writeHeader() error {
	var h [HeaderLen]byte
	copy(h[0:4], headerIdent)
	binary.LittleEndian.PutUint16(h[4:6], headerMagic)
	h[6] = HeaderVersion
	_, err := w.buf.Write(h[:])
	return err
}

// Offset 返回下一条记录写入的起始字节偏移（等于当前文件末尾，含块尾填充）。
// 它是 Reader.OpenOffset 唯一接受的起点，也是上层登记记录边界的依据。
//
// 注意：经 AsyncWriter 异步写入时，Offset 只反映已由写协程消费的部分，
// 调用方不能据此给刚入队的记录记账；那种场景应使用 AsyncWriter.Flushed
// 发布的进度，或事后用 Reader 扫描得到记录边界。
func (w *Writer) Offset() int64 { return w.index.Load() }

// Write 追加一条完整记录。
func (w *Writer) Write(payload []byte) error {
	w.mu.Lock()
	defer w.mu.Unlock()
	if w.closed {
		return ErrClosed
	}
	return w.write(payload)
}

// write 按 LevelDB log 规范把一个逻辑记录切成若干 chunk，调用方需持有 mu。
func (w *Writer) write(p []byte) error {
	space := int(BlockLen - w.Offset()%BlockLen)
	// 块尾放不下一个 chunk 头：补零对齐到下一块
	if space < HeaderLen {
		if err := w.writeZeros(space); err != nil {
			return err
		}
		space = BlockLen
	}
	if len(p)+HeaderLen <= space {
		return w.writeChunk(p, chunkFull)
	}
	// 跨块：FIRST 填满当前块剩余空间，此后每个 MIDDLE 独占整块，最后一段为 LAST
	room := space - HeaderLen
	if err := w.writeChunk(p[:room], chunkFirst); err != nil {
		return err
	}
	used, left := room, len(p)-room
	for left > dataLen {
		if err := w.writeChunk(p[used:used+dataLen], chunkMiddle); err != nil {
			return err
		}
		used += dataLen
		left -= dataLen
	}
	return w.writeChunk(p[used:], chunkLast)
}

// writeChunk 写入一个 chunk 头与 payload，并推进索引。
func (w *Writer) writeChunk(p []byte, typ byte) error {
	var h [HeaderLen]byte
	binary.LittleEndian.PutUint32(h[0:4], crc32.Update(crcSeed[typ], crc32.IEEETable, p))
	binary.LittleEndian.PutUint16(h[4:6], uint16(len(p)))
	h[6] = typ
	if _, err := w.buf.Write(h[:]); err != nil {
		return err
	}
	if len(p) > 0 {
		if _, err := w.buf.Write(p); err != nil {
			return err
		}
	}
	w.index.Add(int64(HeaderLen + len(p)))
	return nil
}

// writeZeros 写入块尾填充；n 一定小于 HeaderLen。
func (w *Writer) writeZeros(n int) error {
	var zeros [HeaderLen]byte
	if _, err := w.buf.Write(zeros[:n]); err != nil {
		return err
	}
	w.index.Add(int64(n))
	return nil
}

// Flush 把用户态缓冲写入操作系统。只有 Flush 成功后，已写入的字节才对
// 并发 Reader 可见——AsyncWriter 依赖这一点发布 (flushedStoreSeq, flushedOffset)。
func (w *Writer) Flush() error {
	w.mu.Lock()
	defer w.mu.Unlock()
	if w.closed {
		return ErrClosed
	}
	return w.buf.Flush()
}

// Close 刷盘并关闭文件。可重复调用，第二次返回 nil。
func (w *Writer) Close() error {
	w.mu.Lock()
	defer w.mu.Unlock()
	if w.closed {
		return nil
	}
	w.closed = true
	err := w.buf.Flush()
	if cerr := w.fp.Close(); err == nil {
		err = cerr
	}
	return err
}

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
		// 与 Python DataStoreReader 的错误语义保持一致：版本不兼容时明确提示升级 SDK。
		return fmt.Errorf("%w: invalid run version: %d (expected %d); "+
			"this file was created by an incompatible version of SwanLab SDK, "+
			"please use a matching SDK version to sync this run", ErrBadHeader, v, HeaderVersion)
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
	r.index += rest
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
