// Package uploader 把 intake 的记录攒批、按类型分轨上传，并提供可等待的排空状态。
//
// 各轨独立重试；列记录先于依赖它的数据记录送达；入口封闭且全部记录终结才判定 Drained。
package uploader

import (
	"context"
	"errors"
	"sync"
	"sync/atomic"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// ErrNotDrained 表示 Close 时仍有记录未终结（配额由调用方决定补传或丢弃）。
var ErrNotDrained = errors.New("uploader: closed before drained")

// ErrClosed 表示 uploader 已停止，不能再接受调度请求。
var ErrClosed = errors.New("uploader: closed")

// Options 是 Uploader 的构造参数；In、Budget、Senders 必填。
type Options struct {
	In        <-chan *recordv1.Record // 必填，通常为 buffer.Intake.UploadOut()
	Budget    *budget.Budget          // 必填，记录入队时已由 intake 预留
	Senders   map[record.Type]Sender  // 必填，至少一条；未注册类型按终结处理并告警
	Batch     int                     // <=0 取 config.DefaultRecordBatch
	Interval  time.Duration           // <=0 取 config.DefaultRecordInterval
	Retry     RetryPolicy             // 零值取 DefaultRetryPolicy()
	Tracker   *Tracker                // 可空，内部创建
	Warn      *WarningThrottle        // 可空，内部创建
	Context   context.Context         // 可空，作为 sender 调用的父 ctx
	OnDrained func()                  // 可空；在判定协程中调用一次，禁止阻塞或再调用 Close
}

// Stats 是上传管线的整体快照。
type Stats struct {
	Buffered        int
	TotalRecords    int64 // Tracker
	UploadedRecords int64 // Tracker
	Tracks          []TrackStats
}

// Uploader 是上传管线门面：独占读入口通道，攒批后分发到类型轨道。
type Uploader struct {
	in       <-chan *recordv1.Record
	batch    int
	interval time.Duration

	budget   *budget.Budget
	tracker  *Tracker
	registry *Registry
	buf      *recordBuffer

	stop        chan struct{}
	finishCh    chan struct{}
	flushReq    chan chan struct{}
	progress    chan struct{}
	done        chan struct{}
	batcherDone chan struct{}

	doneOnce  sync.Once
	closeOnce sync.Once
	wg        sync.WaitGroup

	ctx    context.Context
	cancel context.CancelFunc

	onDrained func()

	buffered    atomic.Int64
	warnedTypes map[record.Type]struct{}
	// warnUnsupported 输出未注册类型的告警；单测可替换。
	warnUnsupported func(typ record.Type, count int)
}

// New 校验参数并启动后台协程；调用方必须在退出前 Close。
func New(opts Options) (*Uploader, error) {
	if opts.In == nil {
		return nil, errors.New("uploader: In is required")
	}
	if opts.Budget == nil {
		return nil, errors.New("uploader: Budget is required")
	}
	if len(opts.Senders) == 0 {
		return nil, errors.New("uploader: at least one Sender is required")
	}

	batch := opts.Batch
	if batch <= 0 {
		batch = config.DefaultRecordBatch
	}
	interval := opts.Interval
	if interval <= 0 {
		interval = time.Duration(config.DefaultRecordInterval * float64(time.Second))
	}
	retry := opts.Retry
	if retry.isZero() {
		retry = DefaultRetryPolicy()
	}
	tracker := opts.Tracker
	if tracker == nil {
		tracker = NewTracker()
	}
	warn := opts.Warn
	if warn == nil {
		warn = NewWarningThrottle()
	}
	parent := opts.Context
	if parent == nil {
		parent = context.Background()
	}
	ctx, cancel := context.WithCancel(parent)

	u := &Uploader{
		in:              opts.In,
		batch:           batch,
		interval:        interval,
		budget:          opts.Budget,
		tracker:         tracker,
		stop:            make(chan struct{}),
		finishCh:        make(chan struct{}, 1),
		flushReq:        make(chan chan struct{}),
		progress:        make(chan struct{}, 1),
		done:            make(chan struct{}),
		batcherDone:     make(chan struct{}),
		ctx:             ctx,
		cancel:          cancel,
		onDrained:       opts.OnDrained,
		warnedTypes:     make(map[record.Type]struct{}),
		warnUnsupported: defaultUnsupportedWarning,
	}
	u.buf = newRecordBuffer(opts.Budget, tracker)
	u.registry = newRegistry(registryConfig{
		senders: opts.Senders,
		batch:   batch,
		retry:   retry,
		budget:  opts.Budget,
		tracker: tracker,
		warn:    warn,
		ctx:     ctx,
		notify:  u.signal,
	})
	u.registry.start()

	u.wg.Add(2)
	go u.runBatcher()
	go u.runDrainWatcher()
	return u, nil
}

// defaultUnsupportedWarning 对未注册类型告警，每种类型只提醒一次。
func defaultUnsupportedWarning(typ record.Type, count int) {
	label := string(typ)
	if label == "" {
		label = "unknown"
	}
	console.Warningf("uploader: no sender registered for %q records, %d record(s) dropped", label, count)
}

// signal 非阻塞地产生一次进度信号，重复信号合并。
func (u *Uploader) signal() {
	select {
	case u.progress <- struct{}{}:
	default:
	}
}

// runBatcher 串行处理入队、定量/定时出批、显式 flush 与收尾。
func (u *Uploader) runBatcher() {
	defer u.wg.Done()
	defer close(u.batcherDone)

	timer := time.NewTimer(u.interval)
	if !timer.Stop() {
		<-timer.C
	}
	armed := false
	disarm := func() {
		if !armed {
			return
		}
		if !timer.Stop() {
			select {
			case <-timer.C:
			default:
			}
		}
		armed = false
	}
	defer disarm()

	for {
		var timerC <-chan time.Time
		if armed {
			timerC = timer.C
		}
		select {
		case rec, ok := <-u.in:
			if !ok {
				disarm()
				u.flush()
				u.registry.markIngestClosed()
				u.signal()
				return
			}
			u.buf.add(rec)
			u.buffered.Store(int64(u.buf.len()))
			if u.buf.len() >= u.batch {
				disarm()
				u.flush()
			} else if !armed {
				timer.Reset(u.interval)
				armed = true
			}
		case <-timerC:
			armed = false
			u.flush()
		case <-u.finishCh:
			disarm()
			u.flush()
		case ack := <-u.flushReq:
			// 先读尽入口通道，保证 flush 覆盖调用前的全部输入。
			closed := u.pumpReady()
			disarm()
			u.flush()
			if closed {
				u.registry.markIngestClosed()
				u.signal()
				close(ack)
				return
			}
			close(ack)
		case <-u.stop:
			return
		}
	}
}

// pumpReady 非阻塞读尽入口通道；返回 true 表示通道已关闭。
func (u *Uploader) pumpReady() (closed bool) {
	for {
		select {
		case rec, ok := <-u.in:
			if !ok {
				return true
			}
			u.buf.add(rec)
			u.buffered.Store(int64(u.buf.len()))
			if u.buf.len() >= u.batch {
				u.flush()
			}
		default:
			return false
		}
	}
}

// flush 取走整个缓冲区并按类型分发。
func (u *Uploader) flush() {
	records := u.buf.drain()
	u.buffered.Store(0)
	if len(records) == 0 {
		return
	}
	u.route(groupRecords(records))
	u.signal()
}

// route 把已注册类型交给轨道，未注册类型就地终结并告警。
func (u *Uploader) route(groups []recordGroup) {
	accepted := make([]recordGroup, 0, len(groups))
	for _, group := range groups {
		if u.registry.has(group.typ) {
			accepted = append(accepted, group)
			continue
		}
		u.terminateUnsupported(group)
	}
	u.registry.dispatch(accepted)
}

// terminateUnsupported 终结未注册类型的记录：归还配额、计数并告警一次。
func (u *Uploader) terminateUnsupported(group recordGroup) {
	for _, rec := range group.records {
		u.budget.Release(1, record.Size(rec))
	}
	u.tracker.AdvanceRecords(int64(len(group.records)))
	if _, warned := u.warnedTypes[group.typ]; warned {
		return
	}
	u.warnedTypes[group.typ] = struct{}{}
	u.warnUnsupported(group.typ, len(group.records))
}

// runDrainWatcher 消费进度信号，排空时关闭 Done 并回调一次。
func (u *Uploader) runDrainWatcher() {
	defer u.wg.Done()
	if u.isDrained() {
		u.markDrained()
		return
	}
	for {
		select {
		case <-u.progress:
			if u.isDrained() {
				u.markDrained()
				return
			}
		case <-u.stop:
			return
		}
	}
}

// isDrained 判定入口封闭、缓冲区为空且所有轨道空闲。
func (u *Uploader) isDrained() bool {
	return u.buffered.Load() == 0 && u.registry.drained()
}

// markDrained 只执行一次：关闭 Done 并回调 OnDrained。
func (u *Uploader) markDrained() {
	u.doneOnce.Do(func() {
		close(u.done)
		if u.onDrained != nil {
			u.onDrained()
		}
	})
}

// RequestFinish 立即出批并打断退避；非阻塞、幂等。
func (u *Uploader) RequestFinish() {
	u.registry.beginFinish()
	select {
	case u.finishCh <- struct{}{}:
	default:
	}
}

// Drained 报告是否已排空。
func (u *Uploader) Drained() bool {
	select {
	case <-u.done:
		return true
	default:
		return false
	}
}

// Done 在排空时关闭。
func (u *Uploader) Done() <-chan struct{} { return u.done }

// Wait 等待排空；超时返回 ctx.Err()，已 Close 返回 ErrNotDrained。
func (u *Uploader) Wait(ctx context.Context) error {
	select {
	case <-u.done:
		return nil
	default:
	}
	select {
	case <-u.done:
		return nil
	case <-u.stop:
		return ErrNotDrained
	case <-ctx.Done():
		return ctx.Err()
	}
}

// Suspend 暂停取批，在途批次照常完成；仅测试使用。
func (u *Uploader) Suspend() { u.registry.suspend() }

// Resume 恢复取批，并把暂停时长顺延到各轨退避截止时间。
func (u *Uploader) Resume() { u.registry.resume() }

// Stats 返回缓冲、记录进度与各轨状态。
func (u *Uploader) Stats() Stats {
	snapshot := u.tracker.Snapshot()
	return Stats{
		Buffered:        int(u.buffered.Load()),
		TotalRecords:    snapshot.TotalRecords,
		UploadedRecords: snapshot.UploadedRecords,
		Tracks:          u.registry.stats(),
	}
}

// PendingRecords 返回仍在轨道队列中、尚未出队的记录。
func (u *Uploader) PendingRecords() []*recordv1.Record { return u.registry.pendingRecords() }

// LastFailures 返回最近一次完成批次的失败集合。
func (u *Uploader) LastFailures() []Failure { return u.registry.lastFailures() }

// Close 强制停止所有协程，不等待队列排空；未排空时返回 ErrNotDrained。
func (u *Uploader) Close() error {
	u.closeOnce.Do(func() {
		close(u.stop)
		u.cancel()
		u.registry.stop()
	})
	u.registry.wait()
	u.wg.Wait()
	if u.Drained() {
		return nil
	}
	return ErrNotDrained
}

// FlushSchedulingForTest 触发一次立即出批并等待分发完成，不等待上传完成。
func (u *Uploader) FlushSchedulingForTest() error {
	ack := make(chan struct{})
	select {
	case u.flushReq <- ack:
	case <-u.batcherDone:
		return ErrClosed
	case <-u.stop:
		return ErrClosed
	}
	select {
	case <-ack:
		return nil
	case <-u.batcherDone:
		return ErrClosed
	case <-u.stop:
		return ErrClosed
	}
}
