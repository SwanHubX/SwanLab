package store

import (
	"bufio"
	"encoding/binary"
	"hash/crc32"
	"os"
	"sync"
	"sync/atomic"
)

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
		if err := w.writeZeros(int(space)); err != nil {
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
