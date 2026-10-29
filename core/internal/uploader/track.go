package uploader

import (
	"context"
	"fmt"
	"sync"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// Result 是一批记录的三态结果：每个输入记录必须恰好出现在 Confirmed 或 Remaining 之一。
//
// Confirmed 表示已终结（含 4xx、House 202、编码失败等确定性失败），Remaining 表示
// 可重试；Failures 供 LastFailures 与确认摘要使用，可跨两类记录。
type Result struct {
	Confirmed []*recordv1.Record
	Remaining []*recordv1.Record
	Failures  []Failure
}

// Failure 描述一次失败：Message 面向用户，Records 是受影响的记录，Err 可为 nil。
type Failure struct {
	Message string
	Records []*recordv1.Record
	Err     error
}

// Sender 处理某一类型的整批记录。实现必须并发安全：同一 Track 串行调用，
// 不同 Track 可能并发。Upload 不得修改入参切片，也不得返回入参之外的记录。
type Sender interface {
	Upload(ctx context.Context, records []*recordv1.Record) Result
}

// TrackStats 是单条轨道的运行状态快照。
type TrackStats struct {
	Type         record.Type
	Queue        int
	InFlight     int
	Confirmed    int64
	FailedRounds int
	RetryAt      time.Time // 零值表示无待重试
}

// track 是同一类型记录的 FIFO 队列与重试状态，字段由 Registry.mu 保护。
type track struct {
	typ       record.Type
	sender    Sender
	queue     []*recordv1.Record
	inflight  int
	confirmed int64
	failed    int
	retryAt   time.Time
	wake      chan struct{}
}

// wakeUp 非阻塞地唤醒轨道 worker，重复信号合并为一次。
func (t *track) wakeUp() {
	select {
	case t.wake <- struct{}{}:
	default:
	}
}

// registryConfig 是 Registry 的构造参数，由 Uploader 填充。
type registryConfig struct {
	senders map[record.Type]Sender
	batch   int
	retry   RetryPolicy
	budget  *budget.Budget
	tracker *Tracker
	warn    *WarningThrottle
	ctx     context.Context
	notify  func()
	now     func() time.Time
}

// Registry 管理所有类型轨道：入队、门控取批、结果记账与停复机。
//
// 全局只有一把 Registry 锁，sender 调用与配额归还都在锁外完成；缓冲区分发、
// 门控判定与出队在同一临界区内，保证被门控的轨道不会越过依赖。
type Registry struct {
	cfg registryConfig

	mu           sync.Mutex
	tracks       map[record.Type]*track
	order        []record.Type
	gates        map[record.Type][]record.Type
	finishing    bool
	ingestClosed bool
	suspendedAt  time.Time
	suspended    bool
	stopped      bool
	failures     []Failure

	stopCh   chan struct{}
	stopOnce sync.Once
	workers  sync.WaitGroup
}

// newRegistry 按固定顺序注册 sender 并建立门控表；未提供 sender 的类型跳过。
func newRegistry(cfg registryConfig) *Registry {
	if cfg.batch <= 0 {
		cfg.batch = config.DefaultRecordBatch
	}
	if cfg.warn == nil {
		cfg.warn = NewWarningThrottle()
	}
	if cfg.ctx == nil {
		cfg.ctx = context.Background()
	}
	if cfg.now == nil {
		cfg.now = time.Now
	}
	r := &Registry{
		cfg:    cfg,
		tracks: make(map[record.Type]*track, len(cfg.senders)),
		gates:  defaultGates(),
		stopCh: make(chan struct{}),
	}
	for _, typ := range dispatchOrder {
		sender, ok := cfg.senders[typ]
		if !ok || sender == nil {
			continue
		}
		r.order = append(r.order, typ)
		r.tracks[typ] = &track{typ: typ, sender: sender, wake: make(chan struct{}, 1)}
	}
	return r
}

// defaultGates 返回固定的跨轨依赖：列是所有数据轨的前置，日志不受门控。
func defaultGates() map[record.Type][]record.Type {
	return map[record.Type][]record.Type{
		record.TypeScalar: {record.TypeColumn},
		record.TypeMedia:  {record.TypeColumn},
		record.TypeSave:   {record.TypeColumn},
	}
}

// start 为每条轨道启动一个 worker。
func (r *Registry) start() {
	for _, typ := range r.order {
		r.workers.Add(1)
		go r.tracks[typ].run(r)
	}
}

// stop 置停止位并关闭停止通道；worker 退出后由 wait 回收。
func (r *Registry) stop() {
	r.stopOnce.Do(func() {
		r.mu.Lock()
		r.stopped = true
		r.mu.Unlock()
		close(r.stopCh)
		for _, typ := range r.order {
			r.tracks[typ].wakeUp()
		}
	})
}

// wait 等待所有 worker 退出。
func (r *Registry) wait() { r.workers.Wait() }

// dispatch 把一批分组记录追加到各自轨道队尾并唤醒；未注册类型忽略。
func (r *Registry) dispatch(groups []recordGroup) {
	if len(groups) == 0 {
		return
	}
	r.mu.Lock()
	awake := make([]*track, 0, len(groups))
	for _, group := range groups {
		t := r.tracks[group.typ]
		if t == nil {
			continue
		}
		t.queue = append(t.queue, group.records...)
		awake = append(awake, t)
	}
	r.mu.Unlock()

	for _, t := range awake {
		t.wakeUp()
	}
	r.notify()
}

// beginFinish 进入收尾：不再等待退避窗口，队列中的记录立即重试。
func (r *Registry) beginFinish() {
	r.mu.Lock()
	r.finishing = true
	for _, typ := range r.order {
		r.tracks[typ].retryAt = time.Time{}
	}
	r.mu.Unlock()
	r.wakeAll()
}

// markIngestClosed 记录入口已封闭；这是 Drained 成立的前提。
func (r *Registry) markIngestClosed() {
	r.mu.Lock()
	r.ingestClosed = true
	r.finishing = true
	for _, typ := range r.order {
		r.tracks[typ].retryAt = time.Time{}
	}
	r.mu.Unlock()
	r.wakeAll()
}

// drained 报告「入口封闭且所有轨道空闲」。调用方还需确认攒批缓冲区为空。
func (r *Registry) drained() bool {
	r.mu.Lock()
	defer r.mu.Unlock()
	if !r.ingestClosed {
		return false
	}
	for _, typ := range r.order {
		t := r.tracks[typ]
		if len(t.queue) > 0 || t.inflight > 0 {
			return false
		}
	}
	return true
}

// suspend 暂停取批：已出队的批次照常完成，队列保留在内存中。
func (r *Registry) suspend() {
	r.mu.Lock()
	if !r.suspended {
		r.suspended = true
		r.suspendedAt = r.cfg.now()
	}
	r.mu.Unlock()
}

// resume 恢复取批，并把暂停时长顺延到各轨的退避截止时间（冻结语义）。
func (r *Registry) resume() {
	r.mu.Lock()
	if !r.suspended {
		r.mu.Unlock()
		return
	}
	frozen := r.cfg.now().Sub(r.suspendedAt)
	for _, typ := range r.order {
		t := r.tracks[typ]
		if !t.retryAt.IsZero() {
			t.retryAt = t.retryAt.Add(frozen)
		}
	}
	r.suspended = false
	r.mu.Unlock()
	r.wakeAll()
}

// stats 返回各轨道的有序快照。
func (r *Registry) stats() []TrackStats {
	r.mu.Lock()
	defer r.mu.Unlock()
	out := make([]TrackStats, 0, len(r.order))
	for _, typ := range r.order {
		t := r.tracks[typ]
		out = append(out, TrackStats{
			Type:         typ,
			Queue:        len(t.queue),
			InFlight:     t.inflight,
			Confirmed:    t.confirmed,
			FailedRounds: t.failed,
			RetryAt:      t.retryAt,
		})
	}
	return out
}

// pendingRecords 返回仍在轨道队列中、尚未出队的记录（不含正在上传的批次）。
func (r *Registry) pendingRecords() []*recordv1.Record {
	r.mu.Lock()
	defer r.mu.Unlock()
	var out []*recordv1.Record
	for _, typ := range r.order {
		out = append(out, r.tracks[typ].queue...)
	}
	return out
}

// lastFailures 返回最近一次完成的批次产生的失败集合。
func (r *Registry) lastFailures() []Failure {
	r.mu.Lock()
	defer r.mu.Unlock()
	if len(r.failures) == 0 {
		return nil
	}
	return append([]Failure(nil), r.failures...)
}

// notify 是进度信号，可为空。
func (r *Registry) notify() {
	if r.cfg.notify != nil {
		r.cfg.notify()
	}
}

// wakeAll 广播唤醒所有轨道；重复信号会被合并。
func (r *Registry) wakeAll() {
	r.mu.Lock()
	tracks := make([]*track, 0, len(r.order))
	for _, typ := range r.order {
		tracks = append(tracks, r.tracks[typ])
	}
	r.mu.Unlock()
	for _, t := range tracks {
		t.wakeUp()
	}
}

// run 是轨道 worker 主循环。
func (t *track) run(r *Registry) {
	defer r.workers.Done()
	for {
		batch, wait, exit := r.take(t)
		if exit {
			return
		}
		if batch != nil {
			r.upload(t, batch)
			continue
		}
		if !t.waitFor(r, wait) {
			return
		}
	}
}

// waitFor 等待下一次唤醒；wait<=0 表示只在有新信号或停止时返回。
func (t *track) waitFor(r *Registry, wait time.Duration) bool {
	if wait <= 0 {
		select {
		case <-t.wake:
			return true
		case <-r.stopCh:
			return false
		}
	}
	timer := time.NewTimer(wait)
	defer timer.Stop()
	select {
	case <-t.wake:
		return true
	case <-timer.C:
		return true
	case <-r.stopCh:
		return false
	}
}

// take 在临界区内判定门控与退避，满足条件时取走队首一批。
//
// 返回 exit 表示 worker 应当退出；batch 非空表示可以上传；否则按 wait
// 等待下一次唤醒（wait<=0 表示只等信号）。
func (r *Registry) take(t *track) (batch []*recordv1.Record, wait time.Duration, exit bool) {
	r.mu.Lock()
	defer r.mu.Unlock()
	if r.stopped {
		return nil, 0, true
	}
	if r.ingestClosed && len(t.queue) == 0 && t.inflight == 0 {
		return nil, 0, true
	}
	if r.suspended || len(t.queue) == 0 || !r.gateOpen(t.typ) {
		return nil, 0, false
	}
	now := r.cfg.now()
	if now.Before(t.retryAt) {
		return nil, t.retryAt.Sub(now), false
	}
	n := min(len(t.queue), r.cfg.batch)
	batch = t.queue[:n:n]
	t.queue = t.queue[n:]
	t.inflight += n
	return batch, 0, false
}

// gateOpen 报告 typ 的依赖轨道是否已全部空闲（队列空且无在途批次）。
func (r *Registry) gateOpen(typ record.Type) bool {
	for _, dep := range r.gates[typ] {
		t := r.tracks[dep]
		if t == nil {
			continue
		}
		if len(t.queue) > 0 || t.inflight > 0 {
			return false
		}
	}
	return true
}

// upload 调用 sender 并记账：confirmed 归还配额并计入进度，remaining 按原顺序写回队头。
func (r *Registry) upload(t *track, batch []*recordv1.Record) {
	result := safeUpload(r.cfg.ctx, t.sender, batch)
	if !conserved(batch, result) {
		result = Result{
			Remaining: batch,
			Failures: []Failure{{
				Message: "sender result does not account for every record",
				Records: batch,
			}},
		}
	}

	// 配额与进度在锁外结算；此时 inflight 尚未归零，Drained 不会被提前判定。
	if len(result.Confirmed) > 0 {
		for _, rec := range result.Confirmed {
			r.cfg.budget.Release(1, record.Size(rec))
		}
		r.cfg.tracker.AdvanceRecords(int64(len(result.Confirmed)))
	}

	r.mu.Lock()
	t.inflight -= len(batch)
	if r.stopped {
		r.mu.Unlock()
		return
	}
	r.failures = result.Failures
	t.confirmed += int64(len(result.Confirmed))
	if len(result.Remaining) == 0 {
		t.failed = 0
		t.retryAt = time.Time{}
	} else {
		t.failed++
		t.queue = append(remainingInBatchOrder(batch, result.Remaining), t.queue...)
		t.retryAt = r.cfg.now().Add(r.cfg.retry.Delay(t.failed))
	}
	retrying := len(result.Remaining) > 0
	r.mu.Unlock()

	if retrying {
		r.cfg.warn.Warn()
	} else {
		r.cfg.warn.Reset()
	}
	// 被门控的轨道可能因本批完成而放行，必须重新求值。
	r.wakeAll()
	r.notify()
}

// safeUpload 兜底 sender 的 panic：整批 remaining，等待下一轮重试。
func safeUpload(ctx context.Context, sender Sender, batch []*recordv1.Record) (result Result) {
	defer func() {
		if recovered := recover(); recovered != nil {
			result = Result{
				Remaining: batch,
				Failures: []Failure{{
					Message: fmt.Sprintf("sender panic: %v", recovered),
					Records: batch,
				}},
			}
		}
	}()
	return sender.Upload(ctx, batch)
}

// conserved 校验结果守恒：每条输入记录恰好出现在 Confirmed 或 Remaining 之一，
// 且没有重复、没有入参之外的记录。
func conserved(batch []*recordv1.Record, result Result) bool {
	if len(result.Confirmed)+len(result.Remaining) != len(batch) {
		return false
	}
	seen := make(map[*recordv1.Record]struct{}, len(batch))
	for _, rec := range result.Confirmed {
		if _, dup := seen[rec]; dup {
			return false
		}
		seen[rec] = struct{}{}
	}
	for _, rec := range result.Remaining {
		if _, dup := seen[rec]; dup {
			return false
		}
		seen[rec] = struct{}{}
	}
	for _, rec := range batch {
		if _, ok := seen[rec]; !ok {
			return false
		}
	}
	return true
}

// remainingInBatchOrder 按原批次顺序过滤出 remaining，避免 sender 返回顺序打乱 FIFO。
func remainingInBatchOrder(batch, remaining []*recordv1.Record) []*recordv1.Record {
	if len(remaining) == 0 {
		return nil
	}
	set := make(map[*recordv1.Record]struct{}, len(remaining))
	for _, rec := range remaining {
		set[rec] = struct{}{}
	}
	kept := make([]*recordv1.Record, 0, len(remaining))
	for _, rec := range batch {
		if _, ok := set[rec]; ok {
			kept = append(kept, rec)
		}
	}
	return kept
}
