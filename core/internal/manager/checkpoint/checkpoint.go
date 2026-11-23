// Package checkpoint 负责周期与停机阶段将连续上传进度及本地落盘水位写入 checkpoint.json。
package checkpoint

import (
	"context"
	"os"
	"sync"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/store"
)

const (
	// defaultInterval 是运行期断点兜底写入间隔，也是慢盘自适应退避的上限。
	defaultInterval = 10 * time.Second
	// defaultMinInterval 是事件通知驱动写断点的最小节流间隔：
	// 普通训练默认 1.5s 一批，2s 节流既能将强杀重复窗口压到 ~1 批，又能避免
	// 高频小批写入时在 NFS 上频繁创建与重命名临时文件。
	defaultMinInterval = 2 * time.Second
	// nfsBackoffFactor 在单次写盘较慢时放大下次允许写入的冷却期，防止压垮慢盘/NFS 元数据面。
	nfsBackoffFactor = 4
)

// Store 抽象断点写入所需的本地落盘水位查询能力（生产环境由 *store.AsyncWriter 实现）。
type Store interface {
	Flushed() (seq uint64, offset int64)
	Index() *store.OffsetIndex
}

// Options 配置 checkpoint 持久化写入器。
type Options struct {
	Ctx            context.Context
	File           string
	RunID          string
	RunFile        string
	Store          Store
	Watermark      func() uint64
	Interval       time.Duration
	MinInterval    time.Duration
	InitialSeq     uint64
	InitialWritten bool
	OnTick         func()
	writeFn        func(path string, cp store.Checkpoint, durable bool) error
}

// Writer 管理 checkpoint.json 的运行期事件节流刷盘与终态同步持久化；nil 接收者上的方法均为安全空操作。
type Writer struct {
	ctx         context.Context
	file        string
	runID       string
	runFile     string
	store       Store
	watermark   func() uint64
	interval    time.Duration
	minInterval time.Duration
	onTick      func()
	write       func(path string, cp store.Checkpoint, durable bool) error

	notify    chan struct{}
	startOnce sync.Once
	stopOnce  sync.Once
	stop      chan struct{}
	done      chan struct{}

	mu         sync.Mutex
	warned     bool
	lastSeq    uint64
	hasWrite   bool
	lastSynced bool
}

func defaultWriteCheckpoint(path string, cp store.Checkpoint, durable bool) error {
	if durable {
		return store.WriteCheckpoint(path, cp)
	}
	return store.WriteCheckpointFast(path, cp)
}

// New 创建断点写入器；当路径、Store、RunID 或 Watermark 缺失时，所有写入操作退化为空操作。
func New(opts Options) *Writer {
	interval := opts.Interval
	if interval <= 0 {
		interval = defaultInterval
	}
	minInterval := opts.MinInterval
	if minInterval <= 0 {
		minInterval = defaultMinInterval
	}
	if minInterval > interval {
		minInterval = interval
	}
	write := opts.writeFn
	if write == nil {
		write = defaultWriteCheckpoint
	}
	return &Writer{
		ctx:         opts.Ctx,
		file:        opts.File,
		runID:       opts.RunID,
		runFile:     opts.RunFile,
		store:       opts.Store,
		watermark:   opts.Watermark,
		interval:    interval,
		minInterval: minInterval,
		onTick:      opts.OnTick,
		write:       write,
		notify:      make(chan struct{}, 1),
		lastSeq:     opts.InitialSeq,
		hasWrite:    opts.InitialWritten,
		lastSynced:  opts.InitialWritten,
	}
}

func (w *Writer) enabled() bool {
	return w != nil && w.file != "" && w.store != nil && w.runID != "" && w.watermark != nil
}

// Notify 在上传确认推进水位后非阻塞唤醒断点写入协程；同一节流窗口内多次通知自动合并。
func (w *Writer) Notify() {
	if !w.enabled() {
		return
	}
	select {
	case w.notify <- struct{}{}:
	default:
	}
}

// Start 启动运行期断点持久化后台协程；未启用时安全跳过。
func (w *Writer) Start() {
	if !w.enabled() {
		return
	}
	w.startOnce.Do(func() {
		var ctxDone <-chan struct{}
		if w.ctx != nil {
			ctxDone = w.ctx.Done()
		}
		w.stop = make(chan struct{})
		w.done = make(chan struct{})
		go w.loop(ctxDone)
	})
}

func (w *Writer) loop(ctxDone <-chan struct{}) {
	var throttleTimer *time.Timer
	var throttleC <-chan time.Time
	var nextAllowed time.Time
	defer func() {
		if throttleTimer != nil {
			throttleTimer.Stop()
		}
		close(w.done)
	}()

	ticker := time.NewTicker(w.interval)
	defer ticker.Stop()

	for {
		select {
		case <-ctxDone:
			return
		case <-w.stop:
			return
		case <-w.notify:
			now := time.Now()
			if now.Before(nextAllowed) {
				if throttleTimer == nil {
					throttleTimer = time.NewTimer(nextAllowed.Sub(now))
					throttleC = throttleTimer.C
				}
				continue
			}
			if throttleTimer != nil {
				throttleTimer.Stop()
				throttleTimer = nil
				throttleC = nil
			}
			nextAllowed = w.runFast(nextAllowed)
		case <-throttleC:
			throttleTimer = nil
			throttleC = nil
			nextAllowed = w.runFast(nextAllowed)
		case <-ticker.C:
			if !time.Now().Before(nextAllowed) {
				nextAllowed = w.runFast(nextAllowed)
			}
		}
	}
}

func (w *Writer) runFast(prev time.Time) time.Time {
	elapsed, attempted := w.persist(false)
	if !attempted {
		return prev
	}
	cooldown := max(w.minInterval, elapsed*nfsBackoffFactor)
	if cooldown > w.interval {
		cooldown = w.interval
	}
	return time.Now().Add(cooldown)
}

// Stop 停止断点协程并等待其退出；未启动或已停止时可重复安全调用。
func (w *Writer) Stop() {
	if w == nil || w.stop == nil {
		return
	}
	w.stopOnce.Do(func() {
		close(w.stop)
	})
	<-w.done
}

// Flush 在停机或排空完成后执行完整 fsync 的最终断点持久化。
func (w *Writer) Flush() {
	w.persist(true)
}

// Persist 把当前已连续确认且已落盘的进度 min(Watermark(), flushedSeq) 同步持久化到 checkpoint.json。
// 未配置 Watermark 时不会写入；写入失败仅告警一次并留待后续重试，不中断会话运行。
func (w *Writer) Persist() {
	w.persist(true)
}

func (w *Writer) persist(durable bool) (time.Duration, bool) {
	if !w.enabled() {
		return 0, false
	}
	w.mu.Lock()
	defer w.mu.Unlock()

	flushedSeq, flushedOff := w.store.Flushed()
	targetSeq := w.watermark()
	if flushedSeq < targetSeq {
		targetSeq = flushedSeq
	}
	if targetSeq == 0 {
		return 0, false
	}
	if w.hasWrite && targetSeq == w.lastSeq && (!durable || w.lastSynced) {
		return 0, false
	}

	var (
		offset int64
		err    error
	)
	if targetSeq == flushedSeq {
		offset = flushedOff
	} else {
		offset, err = store.EndOffset(w.runFile, w.store.Index(), targetSeq, flushedOff)
		if err != nil {
			w.warn(err)
			return 0, false
		}
	}

	cp := store.Checkpoint{
		Version:   store.CheckpointVersion,
		RunID:     w.runID,
		StoreSeq:  targetSeq,
		Offset:    offset,
		UpdatedAt: time.Now(),
	}
	start := time.Now()
	err = w.write(w.file, cp, durable)
	elapsed := time.Since(start)
	if err != nil {
		w.warn(err)
		return elapsed, true
	}
	w.warned = false
	w.lastSeq = targetSeq
	w.hasWrite = true
	w.lastSynced = durable
	if w.onTick != nil {
		w.onTick()
	}
	return elapsed, true
}

// RemoveFile 清理回滚期间产生的断点文件（忽略删除错误）。
func (w *Writer) RemoveFile() {
	if w == nil || w.file == "" {
		return
	}
	_ = os.Remove(w.file)
}

func (w *Writer) warn(err error) {
	if w.warned {
		return
	}
	w.warned = true
	console.Warningf("checkpoint: write %s: %v", w.file, err)
}
