package store

import (
	"context"
	"errors"
	"sync"
	"sync/atomic"
	"time"
)

const (
	defaultFlushBytes    = 64 << 10
	defaultFlushInterval = time.Second
)

// Request 是一条待落盘记录。StoreSeq 由 intake 分配，会话内从 1 开始严格递增，
// 与 Record.num 无关（num 可重复、可覆盖，不能作为进度游标）。
// Payload 的所有权随发送转移给 AsyncWriter，调用方不得复用。
type Request struct {
	StoreSeq uint64
	Payload  []byte
}

// AsyncOptions 控制 AsyncWriter 的入口容量与刷盘节奏。
type AsyncOptions struct {
	// Queue 是写入入口容量（条），必须为正，是背压的唯一来源。
	Queue int
	// FlushBytes 是未刷盘字节数阈值，<=0 取 64 KiB。
	FlushBytes int
	// FlushInterval 是定时刷盘间隔，<=0 取 1s。
	FlushInterval time.Duration
	// OnError 在首次写入或刷盘失败时回调一次，供上层置 storeDegraded；
	// 之后丢弃写入但继续消费入口，避免阻断上传主链路。
	OnError func(error)
}

// flushMark 是已刷盘快照：StoreSeq 及之前的记录已完整落盘，Offset 为其结束偏移。
type flushMark struct {
	StoreSeq uint64
	Offset   int64
}

// AsyncWriter 在独立协程消费写入请求，按定量 / 定时 / Close 三条路径刷盘。
//
// 只有整条记录的全部 chunk 写入且 Flush 成功后，才原子发布 (StoreSeq, Offset)，
// 因此并发 Reader 以此为界不会读到半条记录；发布可能滞后，但不会超前。
type AsyncWriter struct {
	w    *Writer
	ch   chan Request
	ctl  chan control
	done chan struct{}

	flushBytes    int
	flushInterval time.Duration
	onError       func(error)

	mark      atomic.Pointer[flushMark]
	closeOnce sync.Once

	mu  sync.Mutex
	err error
}

// control 是测试用同步请求：flush=false 时只写完已入队记录、不刷盘。
type control struct {
	flush bool
	reply chan error
}

// NewAsyncWriter 启动写入协程；w 的文件头必须已落盘（Create 已保证）。
func NewAsyncWriter(w *Writer, opts AsyncOptions) (*AsyncWriter, error) {
	if w == nil {
		return nil, errors.New("store: nil writer")
	}
	if opts.Queue <= 0 {
		return nil, errors.New("store: queue size must be positive")
	}
	if opts.FlushBytes <= 0 {
		opts.FlushBytes = defaultFlushBytes
	}
	if opts.FlushInterval <= 0 {
		opts.FlushInterval = defaultFlushInterval
	}
	a := &AsyncWriter{
		w:             w,
		ch:            make(chan Request, opts.Queue),
		ctl:           make(chan control),
		done:          make(chan struct{}),
		flushBytes:    opts.FlushBytes,
		flushInterval: opts.FlushInterval,
		onError:       opts.OnError,
	}
	a.mark.Store(&flushMark{Offset: w.Offset()})
	go a.run()
	return a, nil
}

// In 返回写入入口。调用方用 select 实现有界阻塞与 ctx 取消；Close 后不得再发送。
func (a *AsyncWriter) In() chan<- Request { return a.ch }

// Flushed 返回最近一次刷盘后的 (序号, 偏移)，两者始终一致；未刷盘时为 (0, HeaderLen)。
func (a *AsyncWriter) Flushed() (uint64, int64) {
	m := a.mark.Load()
	return m.StoreSeq, m.Offset
}

// Err 返回首次写入或刷盘错误；无错误时为 nil。
func (a *AsyncWriter) Err() error {
	a.mu.Lock()
	defer a.mu.Unlock()
	return a.err
}

// Close 关闭入口，等待排空、刷盘并关闭文件；可重复调用。
func (a *AsyncWriter) Close() error { return a.CloseContext(context.Background()) }

// CloseContext 同 Close，但等待受 ctx 限制；超时后后台继续排空，可再次调用。
func (a *AsyncWriter) CloseContext(ctx context.Context) error {
	a.closeOnce.Do(func() { close(a.ch) })
	select {
	case <-a.done:
		return a.Err()
	case <-ctx.Done():
		return ctx.Err()
	}
}

// FlushForTest 写完已入队记录并立即刷盘，仅供测试。
func (a *AsyncWriter) FlushForTest() error { return a.sync(true) }

// sync 让写协程先排空入口，再按 flush 决定是否刷盘。
func (a *AsyncWriter) sync(flush bool) error {
	reply := make(chan error, 1)
	select {
	case a.ctl <- control{flush: flush, reply: reply}:
	case <-a.done:
		return ErrClosed
	}
	return <-reply
}

func (a *AsyncWriter) run() {
	defer close(a.done)
	ticker := time.NewTicker(a.flushInterval)
	defer ticker.Stop()

	var seq uint64 // 已写入文件的最后一条记录序号
	lastFlush := a.w.Offset()
	for {
		select {
		case req, ok := <-a.ch:
			if !ok {
				a.shutdown(seq)
				return
			}
			seq = a.append(req, seq)
		case <-ticker.C:
			if a.w.Offset() > lastFlush {
				if a.flush(seq) {
					lastFlush = a.w.Offset()
				}
			}
		case c := <-a.ctl:
			var open bool
			seq, open = a.drainAll(seq)
			if !open {
				a.shutdown(seq)
				c.reply <- nil
				return
			}
			var err error
			if c.flush {
				if a.flush(seq) {
					lastFlush = a.w.Offset()
				} else {
					err = a.Err()
				}
			}
			c.reply <- err
		}
		if a.w.Offset()-lastFlush >= int64(a.flushBytes) {
			if a.flush(seq) {
				lastFlush = a.w.Offset()
			}
		}
	}
}

// append 写入一条记录；降级后丢弃，避免阻断上传主链路。
func (a *AsyncWriter) append(req Request, seq uint64) uint64 {
	if a.Err() != nil {
		return seq
	}
	if err := a.w.Write(req.Payload); err != nil {
		a.fail(err)
		return seq
	}
	return req.StoreSeq
}

// drainAll 非阻塞写完已入队记录；入口已关闭时返回 false。
func (a *AsyncWriter) drainAll(seq uint64) (uint64, bool) {
	for {
		select {
		case req, ok := <-a.ch:
			if !ok {
				return seq, false
			}
			seq = a.append(req, seq)
		default:
			return seq, true
		}
	}
}

// flush 刷盘并发布进度。
func (a *AsyncWriter) flush(seq uint64) bool {
	if err := a.w.Flush(); err != nil {
		a.fail(err)
		return false
	}
	a.mark.Store(&flushMark{StoreSeq: seq, Offset: a.w.Offset()})
	return true
}

// shutdown 在入口排空后刷盘、发布最终进度并关闭文件。
func (a *AsyncWriter) shutdown(seq uint64) {
	a.flush(seq)
	if err := a.w.Close(); err != nil {
		a.fail(err)
	}
}

// fail 记录首次错误并回调一次 OnError。
func (a *AsyncWriter) fail(err error) {
	a.mu.Lock()
	first := a.err == nil
	if first {
		a.err = err
	}
	a.mu.Unlock()
	if first && a.onError != nil {
		a.onError(err)
	}
}
