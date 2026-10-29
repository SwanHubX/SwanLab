// Package buffer 实现会话入口分发器 intake：两阶段接收、容量预留与双通道路由。
//
// Reserve 只做只读预检与容量预留；Commit 分配 storeSeq、写本地并分流到上传。
// 本地 run-*.swanlab 是唯一事实源：
//   - 数据类记录同时进入本地与上传；
//   - Start(-1)/Finish(-2) 只落本地；
//   - 崩溃 error_log 只上传；
//   - SAVE_POLICY_END 先落本地、Close 时再上传。
//
// 本地通道满时只做有界反压，调用方保留数据后重试；只有真实写失败才降级为
// 上传兜底。暂停投递期间保留未确认落盘的记录，写盘失败时重新投入上传。
// skip_store 是启动期配置，运行期磁盘写满只置 storeDegraded，不会修改配置。
package buffer

import (
	"context"
	"errors"
	"fmt"
	"sync"
	"sync/atomic"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

const (
	// DefaultQueueWait 是本地通道、上传预算或上传槽位不足时的有界等待上限。
	DefaultQueueWait = 2 * time.Second
	// DefaultUploadQueue 是上传通道的默认容量，与单批上限对齐。
	DefaultUploadQueue = 10_000
	// defaultStorePollInterval 是等待本地或上传通道腾出空间时的轮询间隔。
	defaultStorePollInterval = 10 * time.Millisecond
)

var (
	// ErrClosed 表示 intake 入口已经封闭。
	ErrClosed = errors.New("buffer: intake is closed")
	// ErrOverloaded 表示容量在等待上限内没有释放，调用方保留数据后可重试。
	ErrOverloaded = errors.New("buffer: capacity exhausted")
	// ErrBatchTooLarge 表示单批超过通道容量，重试不会成功。
	ErrBatchTooLarge = errors.New("buffer: batch exceeds queue capacity")
	// ErrOverReserved 表示本次提交超出 Reserve 预留的配额，属于调用方错误。
	ErrOverReserved = errors.New("buffer: commit exceeds reservation")
	// errQueueWait 标记本地通道等待超时；与调用方 ctx 取消区分开。
	errQueueWait = errors.New("buffer: local queue wait exceeded")
)

// Options 是 Intake 的构造参数。
type Options struct {
	// SkipStore 是启动期配置，创建后只读；为 true 时不写本地，Store 必须为 nil。
	SkipStore bool
	// Store 是 AsyncWriter 的写入入口，容量应不小于单批最大记录数。
	Store chan<- store.Request
	// UploadQueue 是上传通道容量，应不小于单批最大记录数；<=0 取 DefaultUploadQueue。
	UploadQueue int
	// Budget 是会话级上传预算。
	Budget *budget.Budget
	// Degraded 是运行期 storeDegraded 标志：写/Flush/序列化失败置位，
	// 容量不足不置位；不会回写配置。为 nil 时内部自建（测试用）。
	Degraded *atomic.Bool
	// Flushed 返回已完整落盘的 storeSeq；!SkipStore 时必填。
	// 用于确认暂停期间的记录已安全落盘，不得回调 Intake。
	Flushed func() uint64
	// QueueWait 是容量不足时的有界等待上限，<=0 取 DefaultQueueWait。
	QueueWait time.Duration
	// StorePollInterval 是本地通道等待的轮询间隔，<=0 取 10ms。
	StorePollInterval time.Duration
}

// Demand 是一次提交对两条通道的容量需求。
type Demand struct {
	// Store 是需要写入本地通道的记录数；SkipStore 时忽略。
	Store int
	// Upload 是需要占用上传预算与通道槽位的记录数（含 SAVE_POLICY_END 记录）。
	Upload int
	// Bytes 是上传记录的序列化字节数估算；Commit 会按实际大小修正，
	// 补扣可能短暂超出上限，调用方应尽量按上界预留。
	Bytes int64
}

// Ticket 是 Reserve 返回的预留凭证，由 Commit 消费或 Release 归还。
type Ticket struct {
	store  int
	upload int
	bytes  int64
	paused bool
}

// pendingSave 是暂存的 SAVE_POLICY_END 记录。
// charged 表示上传配额是否已占用：暂停投递期间暂存的记录在 Close 时才补扣。
type pendingSave struct {
	rec     *recordv1.Record
	charged bool
}

// Intake 是单个会话的入口分发器，并发安全。
type Intake struct {
	skipStore bool
	storeOut  chan<- store.Request
	upload    chan *recordv1.Record
	budget    *budget.Budget
	degraded  *atomic.Bool
	flushed   func() uint64
	wait      time.Duration
	poll      time.Duration
	uploadCap int

	mu            sync.Mutex
	closeMu       sync.Mutex
	closeOnce     sync.Once
	nextSeq       uint64
	pendingStore  int
	pendingUpload int
	pending       []pendingSave
	unconfirmed   []unconfirmedRec
	recovery      []*recordv1.Record
	paused        bool
	pauseSeq      uint64
	closed        bool
}

// New 校验参数并创建 intake。
func New(opts Options) (*Intake, error) {
	if opts.Budget == nil {
		return nil, errors.New("buffer: nil budget")
	}
	if !opts.SkipStore && opts.Store == nil {
		return nil, errors.New("buffer: nil store channel")
	}
	if !opts.SkipStore && opts.Flushed == nil {
		return nil, errors.New("buffer: nil flushed callback")
	}
	if opts.QueueWait <= 0 {
		opts.QueueWait = DefaultQueueWait
	}
	if opts.StorePollInterval <= 0 {
		opts.StorePollInterval = defaultStorePollInterval
	}
	degraded := opts.Degraded
	if degraded == nil {
		degraded = &atomic.Bool{}
	}
	queue := opts.UploadQueue
	if queue <= 0 {
		queue = DefaultUploadQueue
	}
	upload := make(chan *recordv1.Record, queue)
	t := &Intake{
		skipStore: opts.SkipStore,
		storeOut:  opts.Store,
		upload:    upload,
		budget:    opts.Budget,
		degraded:  degraded,
		flushed:   opts.Flushed,
		wait:      opts.QueueWait,
		poll:      opts.StorePollInterval,
		uploadCap: cap(upload),
	}
	return t, nil
}

// UploadOut 返回上传通道；intake 在 Close 时封闭它，消费方读尽缓冲区后退出。
func (t *Intake) UploadOut() <-chan *recordv1.Record { return t.upload }

// Reserve 预留本地通道与上传预算，不修改记录状态。
//
// 预留失败不产生任何副作用。本地通道满时在 QueueWait 内有界等待，仍失败
// 返回 ErrOverloaded，由调用方保留数据后重试。
//
// 本地副本完整但上传预算不足时，本批标记为暂停投递、只落本地；没有完整
// 本地副本时则阻塞等待上传预算，失败返回 ErrOverloaded。
func (t *Intake) Reserve(ctx context.Context, d Demand) (Ticket, error) {
	if d.Store < 0 || d.Upload < 0 || d.Bytes < 0 {
		return Ticket{}, fmt.Errorf("buffer: invalid demand store=%d upload=%d bytes=%d", d.Store, d.Upload, d.Bytes)
	}

	t.mu.Lock()
	if t.closed {
		t.mu.Unlock()
		return Ticket{}, ErrClosed
	}
	paused := t.paused
	if t.degraded.Load() {
		t.collectUnconfirmedLocked()
	}
	t.trimUnconfirmedLocked()
	t.mu.Unlock()

	local := !t.skipStore && !t.degraded.Load()
	tk, err := t.reserveStoreSide(ctx, d, local)
	if err != nil {
		return Ticket{}, err
	}
	if d.Upload == 0 {
		return tk, nil
	}
	return t.reserveUploadSide(ctx, d, local, paused, tk)
}

// reserveStoreSide 预留本地通道；降级或 skip_store 时归还本地预留。
func (t *Intake) reserveStoreSide(ctx context.Context, d Demand, local bool) (Ticket, error) {
	var tk Ticket
	if local && d.Store > 0 {
		err := t.reserveStore(ctx, d.Store)
		switch {
		case err == nil:
			tk.store = d.Store
		case errors.Is(err, ErrBatchTooLarge):
			return Ticket{}, err
		case errors.Is(err, errQueueWait):
			// 本地是唯一事实源：容量不足时反压调用方重试，绝不丢弃/降级。
			return Ticket{}, overload(err)
		case errors.Is(err, context.Canceled), errors.Is(err, context.DeadlineExceeded):
			return Ticket{}, err
		default:
			return Ticket{}, overload(err)
		}
	}
	if !local {
		t.releaseStore(tk.store)
		tk.store = 0
	}
	return tk, nil
}

// reserveUploadSide 预留上传预算与通道槽位；本地副本完整且预算满时转为暂停投递。
func (t *Intake) reserveUploadSide(ctx context.Context, d Demand, local, paused bool, tk Ticket) (Ticket, error) {
	// 暂停投递只在本地副本可靠时成立；降级后必须继续占用上传预算，
	// 否则所有批次都能绕过限额。
	if local && paused {
		tk.paused = true
		return tk, nil
	}
	if !local {
		// 无完整本地副本：禁止暂停投递，只能有界阻塞反压。
		actx, cancel := context.WithTimeoutCause(ctx, t.wait, errQueueWait)
		err := t.budget.Acquire(actx, int64(d.Upload), d.Bytes)
		cause := context.Cause(actx)
		cancel()
		if err != nil {
			t.releaseStore(tk.store)
			if errors.Is(cause, errQueueWait) {
				return Ticket{}, overload(err)
			}
			return Ticket{}, err
		}
	} else if !t.budget.TryAcquire(int64(d.Upload), d.Bytes) {
		t.mu.Lock()
		alreadyPaused := t.paused
		t.paused = true
		t.mu.Unlock()
		if !alreadyPaused {
			console.Warning("upload budget is full, pausing upload delivery; records will stay in the local store")
		}
		tk.paused = true
		return tk, nil
	}

	// 上传通道同样需要预留，保证 Commit 持锁发送时不会阻塞。
	if err := t.reserveUpload(ctx, d.Upload); err != nil {
		t.budget.Release(int64(d.Upload), d.Bytes)
		t.releaseStore(tk.store)
		if errors.Is(err, ErrBatchTooLarge) ||
			errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
			return Ticket{}, err
		}
		return Ticket{}, overload(err)
	}
	tk.upload, tk.bytes = d.Upload, d.Bytes
	return tk, nil
}

// Commit 分配 storeSeq、写本地并分流到上传。票据在返回后即被消费，
// 未使用的预留会自动归还。
//
// 本地已降级时跳过本地写入、继续上传；暂停投递的批次若本地副本已不可用，
// 会改走上传（可能短暂超出预算）以免丢记录。
func (t *Intake) Commit(tk Ticket, records []*recordv1.Record) error {
	if len(records) == 0 {
		t.Release(tk)
		return nil
	}
	for _, rec := range records {
		if rec == nil {
			t.Release(tk)
			return errors.New("buffer: nil record")
		}
	}

	immediate, deferred := routeRecords(records)
	uploadable := len(immediate) + len(deferred)
	payloads, wantStore := t.encodeStore(records)

	t.mu.Lock()
	if t.closed {
		t.mu.Unlock()
		t.Release(tk)
		return ErrClosed
	}
	if wantStore && len(records) > tk.store {
		t.mu.Unlock()
		t.Release(tk)
		return fmt.Errorf("%w: store %d > reserved %d", ErrOverReserved, len(records), tk.store)
	}
	// 票据自身的暂停标记在 ResumeUpload 之后仍然生效。
	paused := t.paused || tk.paused
	if !paused && uploadable > tk.upload {
		t.mu.Unlock()
		t.Release(tk)
		return fmt.Errorf("%w: upload %d > reserved %d", ErrOverReserved, uploadable, tk.upload)
	}
	if !wantStore && tk.store > 0 {
		t.releaseStoreLocked(tk.store)
		tk.store = 0
	}
	t.trimUnconfirmedLocked()

	// 本地通道：预留空间在握，发送不会阻塞。
	seqOf := make([]uint64, len(records))
	t.sendStoreLocked(records, payloads, wantStore, seqOf)
	if tk.store > len(records) {
		t.releaseStoreLocked(tk.store - len(records))
	}

	storeOK := wantStore && !t.degraded.Load()
	if !storeOK && !t.skipStore {
		// 本地已不可用：把此前暂停、尚未确认落盘的记录转入回收队列。
		t.collectUnconfirmedLocked()
	}
	switch {
	case paused && storeOK:
		t.pauseLocked(records, immediate, deferred, tk, seqOf)
	case tk.upload >= len(immediate):
		t.uploadLocked(records, immediate, deferred, tk)
	default:
		// 暂停票据没有上传预留；本地副本又刚失效，只能按条兜底。
		t.fallbackUploadLocked(records, immediate, deferred, tk)
	}
	t.mu.Unlock()

	if !storeOK && !t.skipStore {
		// 本地已不可用，顺带把回收队列推进到上传通道；无空闲槽位时立即返回。
		_, _ = t.DrainRecovery(context.Background())
	}
	return nil
}

// routeRecords 按记录类型拆出立即上传与延迟上传（SAVE_POLICY_END）的下标。
func routeRecords(records []*recordv1.Record) (immediate, deferred []int) {
	immediate = make([]int, 0, len(records))
	deferred = make([]int, 0, len(records))
	for i, rec := range records {
		switch record.TypeOf(rec) {
		case record.TypeStart, record.TypeFinish:
			// 仅落本地。
		case record.TypeSave:
			if rec.GetSave().GetPolicy() == savev1.SavePolicy_SAVE_POLICY_END {
				deferred = append(deferred, i)
				continue
			}
			immediate = append(immediate, i)
		default:
			immediate = append(immediate, i)
		}
	}
	return immediate, deferred
}

// encodeStore 预先序列化落盘 payload；失败则降级并继续上传。
func (t *Intake) encodeStore(records []*recordv1.Record) ([][]byte, bool) {
	if t.skipStore || t.degraded.Load() {
		return nil, false
	}
	payloads := make([][]byte, len(records))
	for i, rec := range records {
		payload, err := proto.Marshal(rec)
		if err != nil {
			t.degrade(fmt.Errorf("marshal record %d for local store: %w", rec.GetNum(), err))
			return nil, false
		}
		payloads[i] = payload
	}
	return payloads, true
}

// sendStoreLocked 把记录写入本地通道并分配 storeSeq；调用方需持有 mu。
func (t *Intake) sendStoreLocked(records []*recordv1.Record, payloads [][]byte, wantStore bool, seqOf []uint64) {
	if !wantStore {
		return
	}
	for i := range records {
		t.nextSeq++
		seqOf[i] = t.nextSeq
		t.storeOut <- store.Request{StoreSeq: t.nextSeq, Payload: payloads[i]}
		t.pendingStore--
	}
}

// uploadLocked 按实际大小修正预算并投入上传通道；需持有 mu。
// 上传槽位已在 Reserve 预留，发送不会阻塞。
func (t *Intake) uploadLocked(records []*recordv1.Record, immediate, deferred []int, tk Ticket) {
	actual := int64(0)
	for _, i := range immediate {
		actual += record.Size(records[i])
	}
	for _, i := range deferred {
		actual += record.Size(records[i])
	}
	t.budget.Adjust(int64(len(immediate)+len(deferred)-tk.upload), actual-tk.bytes)
	for _, i := range immediate {
		t.upload <- records[i]
	}
	for _, i := range deferred {
		t.pending = append(t.pending, pendingSave{rec: records[i], charged: true})
	}
	// 已发出的记录由通道长度计入容量，预留计数整体归还。
	t.releaseUploadLocked(tk.upload)
}

// CommitUploadOnly 提交只上传、不落盘的记录（崩溃 error_log）。
// 没有本地副本，因此在预算内阻塞等待，失败返回错误由调用方重试。
func (t *Intake) CommitUploadOnly(ctx context.Context, rec *recordv1.Record) error {
	if rec == nil {
		return errors.New("buffer: nil record")
	}
	size := record.Size(rec)
	if err := t.budget.Acquire(ctx, 1, size); err != nil {
		if errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
			return err
		}
		return overload(err)
	}
	if err := t.reserveUpload(ctx, 1); err != nil {
		t.budget.Release(1, size)
		if errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
			return err
		}
		return overload(err)
	}

	t.mu.Lock()
	defer t.mu.Unlock()
	if t.closed {
		t.budget.Release(1, size)
		t.releaseUploadLocked(1)
		return ErrClosed
	}
	t.upload <- rec
	t.releaseUploadLocked(1)
	return nil
}

// Release 归还一张未被 Commit 消费的票据；可安全用于零值票据。
func (t *Intake) Release(tk Ticket) {
	if tk.store == 0 && tk.upload == 0 {
		return
	}
	t.mu.Lock()
	defer t.mu.Unlock()
	t.releaseStoreLocked(tk.store)
	if tk.upload > 0 {
		t.budget.Release(int64(tk.upload), tk.bytes)
		t.releaseUploadLocked(tk.upload)
	}
}

// Close 先补投回收队列与暂存的 SAVE_POLICY_END 记录，再关闭上传通道；可重复调用。
// 调用方应在本方法之前关闭本地 store，并在其失败时先调用 DrainRecovery。
func (t *Intake) Close(ctx context.Context) error {
	t.closeMu.Lock()
	defer t.closeMu.Unlock()

	t.mu.Lock()
	t.closed = true
	if t.degraded.Load() {
		t.collectUnconfirmedLocked()
	} else {
		t.unconfirmed = nil
	}
	pending := t.pending
	t.pending = nil
	recovery := t.recovery
	t.recovery = nil
	t.mu.Unlock()

	for i := range recovery {
		size := record.Size(recovery[i])
		t.budget.Adjust(1, size)
		select {
		case t.upload <- recovery[i]:
		case <-ctx.Done():
			t.budget.Adjust(-1, -size)
			t.restoreRecovery(recovery[i:])
			return ctx.Err()
		}
	}
	for i := range pending {
		p := &pending[i]
		if !p.charged {
			if err := t.budget.Acquire(ctx, 1, record.Size(p.rec)); err != nil {
				t.restorePending(pending[i:])
				if errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
					return err
				}
				return overload(err)
			}
			p.charged = true
		}
		select {
		case t.upload <- p.rec:
		case <-ctx.Done():
			t.restorePending(pending[i:])
			return ctx.Err()
		}
	}

	t.closeOnce.Do(func() { close(t.upload) })
	return nil
}

// reserveStore 在有界等待内为 n 条记录预留本地通道空间。
// 容量不足只返回可重试错误；队列满时按 poll 间隔重试。
func (t *Intake) reserveStore(ctx context.Context, n int) error {
	if n > cap(t.storeOut) {
		return fmt.Errorf("%w: batch %d > capacity %d", ErrBatchTooLarge, n, cap(t.storeOut))
	}
	return t.reserveSlots(ctx, n, t.tryReserveStore)
}

// reserveUpload 在有界等待内为 n 条记录预留上传通道槽位。
func (t *Intake) reserveUpload(ctx context.Context, n int) error {
	if n <= 0 {
		return nil
	}
	if n > t.uploadCap {
		return fmt.Errorf("%w: upload batch %d > capacity %d", ErrBatchTooLarge, n, t.uploadCap)
	}
	return t.reserveSlots(ctx, n, t.tryReserveUpload)
}

// reserveSlots 在有界等待内反复尝试占用 n 个槽位。
func (t *Intake) reserveSlots(ctx context.Context, n int, try func(int) bool) error {
	if try(n) {
		return nil
	}
	actx, cancel := context.WithTimeoutCause(ctx, t.wait, errQueueWait)
	defer cancel()
	ticker := time.NewTicker(t.poll)
	defer ticker.Stop()
	for {
		select {
		case <-actx.Done():
			return context.Cause(actx)
		case <-ticker.C:
		}
		if try(n) {
			return nil
		}
	}
}

// tryReserveStore 在空间足够时原子占用 n 个位置。
func (t *Intake) tryReserveStore(n int) bool {
	t.mu.Lock()
	defer t.mu.Unlock()
	if cap(t.storeOut)-len(t.storeOut)-t.pendingStore < n {
		return false
	}
	t.pendingStore += n
	return true
}

// releaseStore 归还本地通道预留；调用方不得持有 mu。
func (t *Intake) releaseStore(n int) {
	if n <= 0 {
		return
	}
	t.mu.Lock()
	defer t.mu.Unlock()
	t.releaseStoreLocked(n)
}

// releaseStoreLocked 归还本地通道预留；调用方需持有 mu。
func (t *Intake) releaseStoreLocked(n int) {
	if n <= 0 {
		return
	}
	if n > t.pendingStore {
		panic("buffer: released more local slots than reserved")
	}
	t.pendingStore -= n
}

func (t *Intake) tryReserveUpload(n int) bool {
	t.mu.Lock()
	defer t.mu.Unlock()
	return t.tryReserveUploadLocked(n)
}

// tryReserveUploadLocked 在槽位足够时原子占用 n 个位置；调用方需持有 mu。
func (t *Intake) tryReserveUploadLocked(n int) bool {
	if n <= 0 {
		return true
	}
	if t.uploadCap-len(t.upload)-t.pendingUpload < n {
		return false
	}
	t.pendingUpload += n
	return true
}

// releaseUploadLocked 归还上传通道槽位；调用方需持有 mu。
func (t *Intake) releaseUploadLocked(n int) {
	if n <= 0 {
		return
	}
	if n > t.pendingUpload {
		panic("buffer: released more upload slots than reserved")
	}
	t.pendingUpload -= n
}

// restorePending 把未发出的延迟记录放回队头，保持原有顺序。
func (t *Intake) restorePending(items []pendingSave) {
	if len(items) == 0 {
		return
	}
	t.mu.Lock()
	t.pending = append(items, t.pending...)
	t.mu.Unlock()
}

// degrade 置为本地存储降级状态并告警一次，只在真实写入失败时调用。
func (t *Intake) degrade(err error) {
	if t.degraded.CompareAndSwap(false, true) {
		console.Error("local store is unavailable, continuing without a local copy:", err)
	}
}

// overload 包装为可重试的容量耗尽错误，保留底层原因。
func overload(err error) error {
	if err == nil {
		return ErrOverloaded
	}
	return fmt.Errorf("%w: %v", ErrOverloaded, err)
}
