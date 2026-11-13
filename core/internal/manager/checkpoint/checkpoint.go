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

const defaultInterval = time.Second

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
	InitialSeq     uint64
	InitialWritten bool
	OnTick         func()
}

// Writer 管理 checkpoint.json 的周期刷盘与终态持久化；nil 接收者上的方法均为安全空操作。
type Writer struct {
	ctx       context.Context
	file      string
	runID     string
	runFile   string
	store     Store
	watermark func() uint64
	interval  time.Duration
	onTick    func()

	startOnce sync.Once
	stopOnce  sync.Once
	stop      chan struct{}
	done      chan struct{}

	mu       sync.Mutex
	warned   bool
	lastSeq  uint64
	hasWrite bool
}

// New 创建断点写入器；当路径、Store、RunID 或 Watermark 缺失时，所有写入操作退化为空操作。
func New(opts Options) *Writer {
	interval := opts.Interval
	if interval <= 0 {
		interval = defaultInterval
	}
	return &Writer{
		ctx:       opts.Ctx,
		file:      opts.File,
		runID:     opts.RunID,
		runFile:   opts.RunFile,
		store:     opts.Store,
		watermark: opts.Watermark,
		interval:  interval,
		onTick:    opts.OnTick,
		lastSeq:   opts.InitialSeq,
		hasWrite:  opts.InitialWritten,
	}
}

func (w *Writer) enabled() bool {
	return w != nil && w.file != "" && w.store != nil && w.runID != "" && w.watermark != nil
}

// Start 启动周期持久化上传断点的后台协程；未启用时安全跳过。
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
		go func() {
			defer close(w.done)
			ticker := time.NewTicker(w.interval)
			defer ticker.Stop()
			for {
				select {
				case <-ctxDone:
					return
				case <-w.stop:
					return
				case <-ticker.C:
					w.Persist()
				}
			}
		}()
	})
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

// Flush 在停机或排空完成后刷写一次最终断点。
func (w *Writer) Flush() {
	w.Persist()
}

// Persist 把当前已连续确认且已落盘的进度 min(Watermark(), flushedSeq) 写入 checkpoint.json。
// 未配置 Watermark 时不会写入；写入失败仅告警一次并留待后续周期重试，不中断会话运行。
func (w *Writer) Persist() {
	if !w.enabled() {
		return
	}
	w.mu.Lock()
	defer w.mu.Unlock()

	flushedSeq, flushedOff := w.store.Flushed()
	targetSeq := w.watermark()
	if flushedSeq < targetSeq {
		targetSeq = flushedSeq
	}
	if targetSeq == 0 {
		return
	}
	if w.hasWrite && targetSeq == w.lastSeq {
		return
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
			return
		}
	}

	cp := store.Checkpoint{
		Version:   store.CheckpointVersion,
		RunID:     w.runID,
		StoreSeq:  targetSeq,
		Offset:    offset,
		UpdatedAt: time.Now(),
	}
	if err = store.WriteCheckpoint(w.file, cp); err != nil {
		w.warn(err)
		return
	}
	w.warned = false
	w.lastSeq = targetSeq
	w.hasWrite = true
	if w.onTick != nil {
		w.onTick()
	}
}

// File 返回断点文件路径。
func (w *Writer) File() string {
	if w == nil {
		return ""
	}
	return w.file
}

// RemoveFile 清理回滚期间产生的断点文件（忽略删除错误）。
func (w *Writer) RemoveFile() {
	if w == nil || w.file == "" {
		return
	}
	_ = os.Remove(w.file)
}

// HasWrite 报告是否已至少成功落盘一次断点（含 Resume 继承的有效断点）。
func (w *Writer) HasWrite() bool {
	if w == nil {
		return false
	}
	w.mu.Lock()
	defer w.mu.Unlock()
	return w.hasWrite
}

// Warned 报告最近一次写入是否处于失败告警状态。
func (w *Writer) Warned() bool {
	if w == nil {
		return false
	}
	w.mu.Lock()
	defer w.mu.Unlock()
	return w.warned
}

func (w *Writer) warn(err error) {
	if w.warned {
		return
	}
	w.warned = true
	console.Warningf("checkpoint: write %s: %v", w.file, err)
}
