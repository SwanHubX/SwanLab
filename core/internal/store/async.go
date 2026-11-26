package store

import (
	"context"
	"errors"
	"sync"
	"sync/atomic"
	"time"
)

// defaultFlushInterval 是定时刷盘间隔：低吞吐时把缓冲里的数据写入。
const defaultFlushInterval = time.Second

// Request 是一条待落盘记录。StoreSeq 由 intake 分配，会话内从 1 开始严格递增，
// 与 Record.num 无关（num 可重复、可覆盖，不能作为进度游标）。
// Payload 的所有权随发送转移给 AsyncWriter，调用方不得复用。
type Request struct {
	StoreSeq int64
	Payload  []byte
}

// AsyncOptions 控制 AsyncWriter 的入口容量与刷盘节奏。
type AsyncOptions struct {
	// Queue 是写入入口容量（条），必须为正，是背压的唯一来源。
	Queue int
	// FlushBytes 是未刷盘字节数阈值；<=0（默认）不按字节阈值刷盘，只靠 bufio 缓冲
	// （32 KiB）写满时自然写出。>0 仅用于测试或需要按量刷盘的调用方。
	FlushBytes int
	// FlushInterval 是定时刷盘间隔，<=0 取 1s。
	FlushInterval time.Duration
	// IndexStride 是 storeSeq→offset 稀疏锚点间隔（条），<=0 取 DefaultIndexStride。
	IndexStride int
	// Index 允许续传时复用 Rebuild 预填充的历史锚点索引；为 nil 时按 IndexStride 新建。
	Index *OffsetIndex
	// ResumeSeq 是续传时已有记录条数，作为初始已刷盘序号；新建会话为 0。
	ResumeSeq int64
	// OnError 在首次写入或刷盘失败时回调一次，供上层置 storeDegraded；
	// 之后丢弃写入但继续消费入口，避免阻断上传主链路。
	OnError func(error)
}

// flushMark 是已交给内核的进度快照：StoreSeq 及之前的记录已完整写出，Offset 为其结束偏移。
type flushMark struct {
	StoreSeq int64
	Offset   int64
}

// AsyncWriter 在独立协程消费写入请求：进度随 bufio 的自然写出推进，显式刷盘只由
// 定时检查、Flush 与 Close 触发（FlushBytes > 0 时再按字节阈值刷盘）。
//
// 发布出去的 (StoreSeq, Offset) 始终对应一条完整记录：刚写完这条记录时，它的末尾
// 可能还在 bufio 缓冲里，此时进度停在上一条已完整写出的记录。进度可能滞后，但不
// 会超前，并发 Reader 以此为界不会读到半条记录。
type AsyncWriter struct {
	w    *Writer
	ch   chan Request
	ctl  chan control
	done chan struct{}

	flushBytes    int
	flushInterval time.Duration
	onError       func(error)
	index         *OffsetIndex

	mark      atomic.Pointer[flushMark]
	closeOnce sync.Once

	mu  sync.Mutex
	err error
}

// control 是写协程的同步控制请求：flush=false 时只写完已入队记录、不刷盘。
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
	if opts.FlushInterval <= 0 {
		opts.FlushInterval = defaultFlushInterval
	}
	idx := opts.Index
	if idx == nil {
		idx = NewOffsetIndex(opts.IndexStride)
	}
	a := &AsyncWriter{
		w:             w,
		ch:            make(chan Request, opts.Queue),
		ctl:           make(chan control),
		done:          make(chan struct{}),
		flushBytes:    opts.FlushBytes,
		flushInterval: opts.FlushInterval,
		onError:       opts.OnError,
		index:         idx,
	}
	a.mark.Store(&flushMark{StoreSeq: opts.ResumeSeq, Offset: w.Offset()})
	go a.run()
	return a, nil
}

// In 返回写入入口。调用方用 select 实现有界阻塞与 ctx 取消；Close 后不得再发送。
func (a *AsyncWriter) In() chan<- Request { return a.ch }

// Flushed 返回已完整交给内核的 (序号, 偏移)，两者严格一致。除显式刷盘外，bufio
// 写满缓冲时的自然写出也会推进进度。新建未写入时为 (0, HeaderLen)，续传时为
// (ResumeSeq, 初始偏移)。
func (a *AsyncWriter) Flushed() (int64, int64) {
	m := a.mark.Load()
	return m.StoreSeq, m.Offset
}

// Index 返回 storeSeq→offset 稀疏索引；写入由写协程串行登记，调用方可并发读取。
func (a *AsyncWriter) Index() *OffsetIndex { return a.index }

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

// Flush 写完已入队记录并刷盘；返回时 Flushed() 覆盖调用前提交的所有请求。
func (a *AsyncWriter) Flush(ctx context.Context) error {
	reply := make(chan error, 1)
	select {
	case a.ctl <- control{flush: true, reply: reply}:
	case <-a.done:
		return ErrClosed
	case <-ctx.Done():
		return ctx.Err()
	}
	select {
	case err := <-reply:
		return err
	case <-ctx.Done():
		return ctx.Err()
	}
}

func (a *AsyncWriter) run() {
	defer close(a.done)
	ticker := time.NewTicker(a.flushInterval)
	defer ticker.Stop()

	mark := a.mark.Load()
	// st 只在写协程内读写，不需要加锁。
	st := writeState{
		seq:         mark.StoreSeq,
		end:         mark.Offset,
		lastFlush:   a.w.Offset(),
		lastWritten: a.w.written(),
	}
	for {
		select {
		case req, ok := <-a.ch:
			if !ok {
				a.shutdown(&st)
				return
			}
			a.append(req, &st)
		case <-ticker.C:
			a.tick(&st)
		case c := <-a.ctl:
			if !a.drainAll(&st) {
				a.shutdown(&st)
				c.reply <- nil
				return
			}
			var err error
			if c.flush {
				if !a.flush(&st) {
					err = a.Err()
				}
			} else {
				// 同步屏障返回前应用字节阈值策略，避免调用方观察到尚未发生的自动刷盘。
				a.autoFlush(&st)
			}
			c.reply <- err
		}
		a.autoFlush(&st)
	}
}

// autoFlush 在未刷盘字节达到阈值时刷盘；失败后进入降级状态。
func (a *AsyncWriter) autoFlush(st *writeState) {
	if st.failed || a.flushBytes <= 0 || a.w.Offset()-st.lastFlush < int64(a.flushBytes) {
		return
	}
	a.flush(st)
}

// writeState 是写协程的私有状态：seq 与 end 是最后一条已写入记录的序号和结束偏移，
// lastFlush 是字节阈值基线，lastWritten 是上次定时检查时已交给内核的字节数。
type writeState struct {
	seq         int64
	end         int64
	lastFlush   int64
	lastWritten int64
	failed      bool // 已降级：继续消费入口但丢弃写入，避免阻断上传主链路
}

// append 写入一条记录并推进进度；降级后只丢弃写入，避免阻断上传主链路。
func (a *AsyncWriter) append(req Request, st *writeState) {
	if st.failed {
		return
	}
	start, end, written, err := a.w.writeRecord(req.Payload)
	if err != nil {
		st.failed = true
		a.fail(err)
		return
	}
	a.index.Record(req.StoreSeq, start, end)
	prevSeq, prevEnd := st.seq, st.end
	st.seq, st.end = req.StoreSeq, end
	// bufio 空间不够时会先写出旧缓冲，所以刚写完的这条记录末尾可能还在缓冲里，
	// 此时进度只能推进到上一条完整写出的记录。
	switch {
	case end <= written:
		a.publish(req.StoreSeq, end)
	case prevSeq > 0 && prevEnd <= written:
		a.publish(prevSeq, prevEnd)
	}
}

// drainAll 非阻塞写完已入队记录；入口已关闭时返回 false。
func (a *AsyncWriter) drainAll(st *writeState) bool {
	for {
		select {
		case req, ok := <-a.ch:
			if !ok {
				return false
			}
			a.append(req, st)
		default:
			return true
		}
	}
}

// tick 是定时刷盘：本周期内 bufio 已自然写出过数据就跳过，避免打断 32 KiB 缓冲
// 产生碎片写入；只有缓冲非空且本周期没有自然写出时才 Flush。
func (a *AsyncWriter) tick(st *writeState) {
	if st.failed {
		return
	}
	written := a.w.written()
	if written > st.lastWritten {
		st.lastWritten = written
		return
	}
	if a.w.Offset() > written {
		a.flush(st)
	}
}

// flush 把缓冲写入内核并发布进度。
func (a *AsyncWriter) flush(st *writeState) bool {
	if err := a.w.Flush(); err != nil {
		st.failed = true
		a.fail(err)
		return false
	}
	off := a.w.Offset()
	a.publish(st.seq, off)
	st.lastFlush = off
	st.lastWritten = off
	return true
}

// publish 推进进度；只在序号或偏移前进时更新，保证 Flushed() 不回退。
func (a *AsyncWriter) publish(seq int64, off int64) {
	old := a.mark.Load()
	if seq < old.StoreSeq || (seq == old.StoreSeq && off <= old.Offset) {
		return
	}
	a.mark.Store(&flushMark{StoreSeq: seq, Offset: off})
}

// shutdown 在入口排空后刷盘、发布最终进度并关闭文件。
func (a *AsyncWriter) shutdown(st *writeState) {
	a.flush(st)
	if err := a.w.Close(); err != nil {
		st.failed = true
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
