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

// Request 是一条待落盘记录。
//
// StoreSeq 是会话内严格单调递增的逻辑写入序号（从 1 开始），由 intake 分配；
// 它与 Record.num 无关：num 会被同号覆盖、-1/-2 只写本地、重复 Save 可多次出现，
// 不能作为落盘进度游标。Payload 是序列化后的 Record，其所有权随发送转移给
// AsyncWriter，调用方不得在发送后复用该切片。
type Request struct {
	StoreSeq uint64
	Payload  []byte
}

// AsyncOptions 控制 AsyncWriter 的入口容量与刷盘节奏。
type AsyncOptions struct {
	// Queue 是写入入口通道容量（条），必须为正——它是有界背压的唯一来源。
	Queue int
	// FlushBytes 是累计未刷盘字节数的刷盘阈值，<=0 时取 64 KiB。
	FlushBytes int
	// FlushInterval 是定时刷盘间隔，<=0 时取 1s。
	FlushInterval time.Duration
	// OnError 在首次写入或刷盘失败时回调一次，供上层把 storeDegraded 置位。
	// 置位后 AsyncWriter 会丢弃后续写入（磁盘副本已不完整），但继续消费入口，
	// 以免阻断上传主链路。
	OnError func(error)
}

// flushMark 是「已刷盘」快照：StoreSeq 及之前的所有记录已完整落盘，
// Offset 是这些记录在文件中的结束偏移。
type flushMark struct {
	StoreSeq uint64
	Offset   int64
}

// AsyncWriter 在独立协程中消费写入请求，按定量/定时/Close 三条路径刷盘。
//
// 进度发布语义：只有一条记录的全部物理分块都写入、且显式 Flush 成功后，
// 才原子发布新的 (StoreSeq, Offset)。因此并发 Reader 以 Offset 为上界时，
// 绝不会读到半条记录；发布可能滞后于实际文件内容，但不会超前。
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

// control 是测试用的同步请求：flush=false 时只把已入队记录写完而不刷盘，
// 用于确定性地验证定量刷盘阈值。
type control struct {
	flush bool
	reply chan error
}

// NewAsyncWriter 启动写入协程。w 的文件头必须已经落盘（Create 已保证）。
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

// In 返回写入入口。调用方通过 select 在该通道上实现有界阻塞与 ctx 取消；
// Close 之后不得再向其发送，否则会 panic on closed channel。
func (a *AsyncWriter) In() chan<- Request { return a.ch }

// Flushed 返回最近一次成功刷盘后的 (逻辑序号, 文件偏移)，两者始终一致。
// 尚未刷出任何记录时返回 (0, HeaderLen)。
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

// Close 关闭写入入口，等待已入队记录全部落盘、刷盘并关闭文件。
// 可重复调用：第二次调用只等待同一排空结果。
func (a *AsyncWriter) Close() error { return a.CloseContext(context.Background()) }

// CloseContext 与 Close 相同，但等待排空受 ctx 限制。
// ctx 到期返回 ctx.Err()，后台继续排空，可再次调用以等待最终结果。
func (a *AsyncWriter) CloseContext(ctx context.Context) error {
	a.closeOnce.Do(func() { close(a.ch) })
	select {
	case <-a.done:
		return a.Err()
	case <-ctx.Done():
		return ctx.Err()
	}
}

// FlushForTest 把已入队记录全部写完并立即刷盘，返回刷盘错误。
// 仅供测试确定性地推进定时刷盘路径，不要在业务代码中调用。
func (a *AsyncWriter) FlushForTest() error { return a.sync(true) }

// sync 请求写入协程处理 ctl：先写完当前已入队记录，再按 flush 决定是否刷盘。
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
	var lastFlush int64 = a.w.Offset()
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

// append 写入一条记录；已降级（存在错误）时丢弃，避免阻断上传主链路。
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

// drainAll 非阻塞地写完当前已入队的全部记录；入口已关闭时返回 false。
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

// flush 刷盘并发布进度，返回是否成功。
func (a *AsyncWriter) flush(seq uint64) bool {
	if err := a.w.Flush(); err != nil {
		a.fail(err)
		return false
	}
	a.mark.Store(&flushMark{StoreSeq: seq, Offset: a.w.Offset()})
	return true
}

// shutdown 在排空入口后收尾：刷盘、发布最终进度并关闭文件。
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
