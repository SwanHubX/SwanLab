// Package uploader 把 intake 的记录攒批并按固定类型顺序上传，提供可等待的排空状态。
// 每个 Session 只有一个上传 worker；统一 pending、退避与告警状态。
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
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
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
	Now       func() time.Time        // 可空，调度时钟，测试可注入
	OnDrained func()                  // 可空；在 worker 中调用一次，禁止阻塞或再调用 Close
}

// Stats 是上传管线的诊断快照；Pending 包含排队与在途记录。
type Stats struct {
	Buffered        int
	Pending         int
	TotalRecords    int64
	UploadedRecords int64
}

type controlRequest struct {
	suspend bool
	ack     chan struct{}
}

// Uploader 是上传管线门面：独占读入口通道，并由单个 worker 完成调度与上传。
type Uploader struct {
	in       <-chan *recordv1.Record
	batch    int
	interval time.Duration
	retry    RetryPolicy
	now      func() time.Time

	budget   *budget.Budget
	tracker  *Tracker
	registry *Registry
	buf      *recordBuffer
	warn     *WarningThrottle

	stop       chan struct{}
	finishCh   chan struct{}
	flushReq   chan chan struct{}
	controlReq chan controlRequest
	done       chan struct{}
	workerDone chan struct{}

	doneOnce  sync.Once
	closeOnce sync.Once
	wg        sync.WaitGroup

	ctx    context.Context
	cancel context.CancelFunc

	onDrained func()

	buffered atomic.Int64
	pending  atomic.Int64

	failuresMu sync.RWMutex
	failures   []Failure

	warnedTypes map[record.Type]struct{}
	// 测试可替换告警出口。
	warnUnsupported func(typ record.Type, count int)
	warnFailure     func(message string)
}

// New 校验参数并启动上传 worker；调用方必须在退出前 Close。
func New(opts Options) (*Uploader, error) {
	if opts.In == nil {
		return nil, errors.New("uploader: In is required")
	}
	if opts.Budget == nil {
		return nil, errors.New("uploader: Budget is required")
	}
	registry := newRegistry(opts.Senders)
	if registry.len() == 0 {
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
	now := opts.Now
	if now == nil {
		now = time.Now
	}
	ctx, cancel := context.WithCancel(parent)

	u := &Uploader{
		in:              opts.In,
		batch:           batch,
		interval:        interval,
		retry:           retry,
		now:             now,
		budget:          opts.Budget,
		tracker:         tracker,
		registry:        registry,
		warn:            warn,
		stop:            make(chan struct{}),
		finishCh:        make(chan struct{}, 1),
		flushReq:        make(chan chan struct{}),
		controlReq:      make(chan controlRequest),
		done:            make(chan struct{}),
		workerDone:      make(chan struct{}),
		ctx:             ctx,
		cancel:          cancel,
		onDrained:       opts.OnDrained,
		warnedTypes:     make(map[record.Type]struct{}),
		warnUnsupported: defaultUnsupportedWarning,
		warnFailure:     func(message string) { console.Warning(message) },
	}
	u.buf = newRecordBuffer(opts.Budget, tracker)
	u.wg.Add(1)
	go u.run()
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

// workerState 只由单个上传 goroutine 持有，统计快照通过 Uploader 的原子计数导出。
type workerState struct {
	u            *Uploader
	pending      []*recordv1.Record
	retryAt      time.Time
	batchAt      time.Time
	suspendedAt  time.Time
	flushAck     chan struct{}
	failed       int
	ingestClosed bool
	finishing    bool
	suspended    bool
}

// run 串行处理入队、出批、上传、统一退避、停复机与排空。
func (u *Uploader) run() {
	defer u.wg.Done()
	defer close(u.workerDone)
	worker := workerState{u: u}
	worker.run()
}

func (w *workerState) run() {
	for {
		if w.drained() {
			w.u.markDrained()
			return
		}
		if w.canDispatch() {
			if w.u.ctx.Err() != nil {
				return
			}
			w.dispatchPending()
			continue
		}
		if w.prepareBatch() {
			continue
		}

		timer, timerC := w.timer()
		exit := false
		select {
		case rec, ok := <-w.input():
			w.handleInput(rec, ok)
		case <-timerC:
			w.handleTimer()
		case <-w.u.finishCh:
			w.handleFinish()
		case ack := <-w.u.flushReq:
			w.handleFlush(ack)
		case req := <-w.u.controlReq:
			w.handleControl(req)
		case <-w.u.stop:
			exit = true
		case <-w.u.ctx.Done():
			exit = true
		}
		if timer != nil {
			timer.Stop()
		}
		if exit {
			return
		}
	}
}

func (w *workerState) drained() bool {
	return w.ingestClosed && w.u.buf.len() == 0 && len(w.pending) == 0
}

func (w *workerState) canDispatch() bool {
	return !w.suspended && len(w.pending) > 0 && w.retryAt.IsZero()
}

func (w *workerState) dispatchPending() {
	var retry bool
	w.pending, retry = w.u.dispatch(w.pending)
	if w.flushAck != nil {
		close(w.flushAck)
		w.flushAck = nil
	}
	if retry {
		w.failed++
		w.retryAt = w.u.now().Add(w.u.retry.Delay(w.failed))
		return
	}
	w.failed = 0
	w.retryAt = time.Time{}
}

func (w *workerState) prepareBatch() bool {
	if w.suspended || len(w.pending) > 0 || w.u.buf.len() == 0 {
		return false
	}
	now := w.u.now()
	if w.batchAt.IsZero() {
		w.batchAt = now.Add(w.u.interval)
	}
	if now.Before(w.batchAt) {
		return false
	}
	w.pending = w.u.drainBuffer()
	w.batchAt = time.Time{}
	return true
}

func (w *workerState) input() <-chan *recordv1.Record {
	if w.ingestClosed || len(w.pending) > 0 {
		return nil
	}
	return w.u.in
}

func (w *workerState) timer() (*time.Timer, <-chan time.Time) {
	if w.suspended {
		return nil, nil
	}
	now := w.u.now()
	switch {
	case len(w.pending) > 0 && !w.retryAt.IsZero():
		timer := time.NewTimer(max(w.retryAt.Sub(now), 0))
		return timer, timer.C
	case len(w.pending) == 0 && w.u.buf.len() > 0 && !w.batchAt.IsZero():
		timer := time.NewTimer(max(w.batchAt.Sub(now), 0))
		return timer, timer.C
	default:
		return nil, nil
	}
}

func (w *workerState) handleInput(rec *recordv1.Record, ok bool) {
	if !ok {
		w.ingestClosed = true
		w.batchAt = time.Time{}
		if w.u.buf.len() > 0 {
			w.pending = w.u.drainBuffer()
		}
		return
	}
	w.u.buf.add(rec)
	w.u.buffered.Store(int64(w.u.buf.len()))
	if w.batchAt.IsZero() {
		w.batchAt = w.u.now().Add(w.u.interval)
	}
	if w.u.buf.len() >= w.u.batch {
		w.pending = w.u.drainBuffer()
		w.batchAt = time.Time{}
	}
}

func (w *workerState) handleTimer() {
	if len(w.pending) > 0 {
		w.retryAt = time.Time{}
		return
	}
	if w.u.buf.len() > 0 {
		w.pending = w.u.drainBuffer()
		w.batchAt = time.Time{}
	}
}

func (w *workerState) handleFinish() {
	if !w.finishing {
		w.finishing = true
		w.retryAt = time.Time{}
	}
	if len(w.pending) == 0 && w.u.buf.len() > 0 {
		w.pending = w.u.drainBuffer()
		w.batchAt = time.Time{}
	}
}

func (w *workerState) handleFlush(ack chan struct{}) {
	if len(w.pending) == 0 && !w.ingestClosed {
		w.ingestClosed = w.u.pumpReady(&w.batchAt)
	}
	if w.suspended || len(w.pending) > 0 || w.u.buf.len() == 0 {
		close(ack)
		return
	}
	w.pending = w.u.drainBuffer()
	w.batchAt = time.Time{}
	w.flushAck = ack
}

func (w *workerState) handleControl(req controlRequest) {
	switch {
	case req.suspend && !w.suspended:
		w.suspended = true
		w.suspendedAt = w.u.now()
	case !req.suspend && w.suspended:
		frozen := w.u.now().Sub(w.suspendedAt)
		if !w.retryAt.IsZero() {
			w.retryAt = w.retryAt.Add(frozen)
		}
		if !w.batchAt.IsZero() {
			w.batchAt = w.batchAt.Add(frozen)
		}
		w.suspended = false
	}
	close(req.ack)
}

// pumpReady 非阻塞读尽当前可读的输入，达到批上限时留给下一轮处理。
func (u *Uploader) pumpReady(batchAt *time.Time) bool {
	for u.buf.len() < u.batch {
		select {
		case rec, ok := <-u.in:
			if !ok {
				return true
			}
			u.buf.add(rec)
			u.buffered.Store(int64(u.buf.len()))
			if batchAt.IsZero() {
				*batchAt = u.now().Add(u.interval)
			}
		default:
			return false
		}
	}
	return false
}

// drainBuffer 把缓冲区移交给 worker 的 pending。
func (u *Uploader) drainBuffer() []*recordv1.Record {
	records := u.buf.drain()
	u.buffered.Store(0)
	u.pending.Store(int64(len(records)))
	return records
}

// dispatch 按固定类型顺序处理一批；失败时返回该组剩余及全部后续组。
func (u *Uploader) dispatch(records []*recordv1.Record) ([]*recordv1.Record, bool) {
	groups := groupRecords(records)
	var failures []Failure
	for i, group := range groups {
		tail := flattenGroups(groups[i+1:])
		u.pending.Store(int64(len(group.records) + len(tail)))
		sender := u.registry.get(group.typ)
		if sender == nil {
			u.terminateUnsupported(group)
			u.pending.Store(int64(len(tail)))
			continue
		}

		result := safeUpload(u.ctx, sender, group.records)
		if !conserved(group.records, result) {
			result = Result{
				Remaining: group.records,
				Failures: []Failure{{
					Message: "sender result does not account for every record",
					Records: group.records,
				}},
			}
		}
		u.confirm(result.Confirmed)
		failures = append(failures, result.Failures...)

		if len(result.Remaining) > 0 {
			pending := remainingInBatchOrder(group.records, result.Remaining)
			pending = append(pending, tail...)
			u.pending.Store(int64(len(pending)))
			u.setFailures(failures)
			u.warn.Warn()
			return pending, true
		}
		if len(result.Failures) > 0 {
			u.warnTerminal(result.Failures)
		}
		u.pending.Store(int64(len(tail)))
	}
	u.pending.Store(0)
	u.setFailures(failures)
	u.warn.Reset()
	return nil, false
}

// flattenGroups 保持固定组序和组内记录顺序。
func flattenGroups(groups []recordGroup) []*recordv1.Record {
	var records []*recordv1.Record
	for _, group := range groups {
		records = append(records, group.records...)
	}
	return records
}

// confirm 终结记录并归还配额。
func (u *Uploader) confirm(records []*recordv1.Record) {
	for _, rec := range records {
		u.budget.Release(1, record.Size(rec))
	}
	u.tracker.AdvanceRecords(int64(len(records)))
}

// terminateUnsupported 终结未注册类型的记录：归还配额、计数并告警一次。
func (u *Uploader) terminateUnsupported(group recordGroup) {
	u.confirm(group.records)
	if _, warned := u.warnedTypes[group.typ]; warned {
		return
	}
	u.warnedTypes[group.typ] = struct{}{}
	u.warnUnsupported(group.typ, len(group.records))
}

// warnTerminal 输出不会被重试的确定性失败，避免静默丢弃数据。
func (u *Uploader) warnTerminal(failures []Failure) {
	for _, failure := range failures {
		u.warnFailure(failure.Message)
	}
}

// setFailures 保存最近一次完整调度尝试的失败集合。
func (u *Uploader) setFailures(failures []Failure) {
	u.failuresMu.Lock()
	u.failures = append(u.failures[:0], failures...)
	u.failuresMu.Unlock()
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

// RequestFinish 立即出批，并仅在首次进入收尾时清除现有退避；非阻塞、幂等。
func (u *Uploader) RequestFinish() {
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

// Wait 等待排空；超时返回 ctx.Err()，worker 未排空退出时返回 ErrNotDrained。
func (u *Uploader) Wait(ctx context.Context) error {
	select {
	case <-u.done:
		return nil
	default:
	}
	select {
	case <-u.done:
		return nil
	case <-u.workerDone:
		if u.Drained() {
			return nil
		}
		return ErrNotDrained
	case <-ctx.Done():
		return ctx.Err()
	}
}

// Suspend 暂停出批与重试，在途调用完成后生效；仅测试使用。
func (u *Uploader) Suspend() { u.control(true) }

// Resume 恢复调度，并把暂停时长顺延到批量与退避截止时间。
func (u *Uploader) Resume() { u.control(false) }

// control 等待 worker 应用停复机状态，避免返回后再产生新上传。
func (u *Uploader) control(suspend bool) {
	req := controlRequest{suspend: suspend, ack: make(chan struct{})}
	select {
	case u.controlReq <- req:
	case <-u.workerDone:
		return
	case <-u.stop:
		return
	}
	select {
	case <-req.ack:
	case <-u.workerDone:
	case <-u.stop:
	}
}

// Stats 返回缓冲、pending 与记录进度。
func (u *Uploader) Stats() Stats {
	snapshot := u.tracker.Snapshot()
	return Stats{
		Buffered:        int(u.buffered.Load()),
		Pending:         int(u.pending.Load()),
		TotalRecords:    snapshot.TotalRecords,
		UploadedRecords: snapshot.UploadedRecords,
	}
}

// Snapshot 返回完整 OperationStats，供 manager 接入 GetOperationStats。
func (u *Uploader) Snapshot() *operationv1.OperationStats { return u.tracker.Snapshot() }

// LastFailures 返回最近一次完整调度尝试的失败集合。
func (u *Uploader) LastFailures() []Failure {
	u.failuresMu.RLock()
	defer u.failuresMu.RUnlock()
	return append([]Failure(nil), u.failures...)
}

// Close 强制停止 worker，不等待队列排空；未排空时返回 ErrNotDrained。
func (u *Uploader) Close() error {
	u.closeOnce.Do(func() {
		close(u.stop)
		u.cancel()
	})
	u.wg.Wait()
	if u.Drained() {
		return nil
	}
	return ErrNotDrained
}

// FlushSchedulingForTest 触发一次立即出批，并等待本轮上传尝试结束。
func (u *Uploader) FlushSchedulingForTest() error {
	ack := make(chan struct{})
	select {
	case u.flushReq <- ack:
	case <-u.workerDone:
		return ErrClosed
	case <-u.stop:
		return ErrClosed
	}
	select {
	case <-ack:
		return nil
	case <-u.workerDone:
		return ErrClosed
	case <-u.stop:
		return ErrClosed
	}
}
