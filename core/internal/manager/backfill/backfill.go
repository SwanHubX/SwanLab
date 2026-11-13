// Package backfill 承载预算暂停后的本地补投：从 pauseSeq 读取已落盘记录，按序投入上传通道，
// 追平后在 intake 锁内恢复实时投递。
package backfill

import (
	"context"
	"errors"
	"fmt"
	"io"
	"sync"
	"sync/atomic"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

const (
	// batchMax 限制单次读取与补投的记录数，避免大 record_batch 占用过多内存。
	batchMax = 1024
	// windowMax 限制交接时持锁读取的尾部记录数。
	windowMax = 256
	// pollInterval 是等待刷盘推进或上传槽位释放的轮询间隔。
	pollInterval = 50 * time.Millisecond
	// retryBase / retryMax 是失败重试的初始与最大退避。
	retryBase = 200 * time.Millisecond
	retryMax  = 2 * time.Second
)

var (
	// errPending 表示本次没有可推进的工作，稍后重试即可。
	errPending = errors.New("backfill: pending")
	// errCorrupt 表示本地记录无法解析，重试不会成功。
	errCorrupt = errors.New("backfill: stopped on corrupt store")
)

// Cursor 指向下一条待读记录：Seq 为已读取的最后一条 storeSeq。
type Cursor struct {
	Seq    uint64
	Offset int64
}

type state struct {
	pauseSeq uint64
	cursor   Cursor
}

// Options 配置本地补投工作器。
type Options struct {
	Ctx           context.Context
	Intake        *buffer.Intake
	Store         *store.AsyncWriter
	RunFile       string
	BatchSize     int32
	Degraded      *atomic.Bool
	IsRunning     func() bool
	IsSuspended   func() bool
	OnStored      func(rec *recordv1.Record, storeSeq uint64)
	InitialPause  uint64
	InitialCursor Cursor
}

// Worker 管理暂停区间的后台补投协程与停机阶段同步补投；nil 接收者上的方法均为安全空操作。
type Worker struct {
	ctx         context.Context
	intake      *buffer.Intake
	store       *store.AsyncWriter
	runFile     string
	batchSize   int
	degraded    *atomic.Bool
	isRunning   func() bool
	isSuspended func() bool
	onStored    func(rec *recordv1.Record, storeSeq uint64)

	wake      chan struct{}
	stop      chan struct{}
	done      chan struct{}
	startOnce sync.Once
	stopOnce  sync.Once
	mu        sync.Mutex
	cancel    context.CancelFunc
	state     state
}

// New 创建补投工作器；当 Intake 或 Store 为 nil 时，所有操作安全退化为空操作。
func New(opts Options) *Worker {
	ctx := opts.Ctx
	if ctx == nil {
		ctx = context.Background()
	}
	batch := int(opts.BatchSize)
	if batch <= 0 {
		batch = batchMax
	} else {
		batch = min(batch, batchMax)
	}
	return &Worker{
		ctx:         ctx,
		intake:      opts.Intake,
		store:       opts.Store,
		runFile:     opts.RunFile,
		batchSize:   batch,
		degraded:    opts.Degraded,
		isRunning:   opts.IsRunning,
		isSuspended: opts.IsSuspended,
		onStored:    opts.OnStored,
		wake:        make(chan struct{}, 1),
		state: state{
			pauseSeq: opts.InitialPause,
			cursor:   opts.InitialCursor,
		},
	}
}

// Start 启动后台补投协程；可重复安全调用，但仅首次生效。
func (w *Worker) Start() {
	if w == nil || w.intake == nil || w.store == nil {
		return
	}
	w.startOnce.Do(func() {
		w.stop = make(chan struct{})
		w.done = make(chan struct{})
		go w.loop()
	})
}

// Stop 取消当前补投尝试并等待协程退出；未启动或已停止时可重复安全调用。
func (w *Worker) Stop() {
	if w == nil || w.stop == nil {
		return
	}
	w.stopOnce.Do(func() {
		close(w.stop)
	})
	w.CancelAttempt()
	<-w.done
}

// Wake 唤醒补投协程；可在任意协程安全调用。
func (w *Worker) Wake() {
	if w == nil || w.wake == nil {
		return
	}
	select {
	case w.wake <- struct{}{}:
	default:
	}
}

// CancelAttempt 取消并清除当前补投尝试；健康状态切换或停止时调用。
func (w *Worker) CancelAttempt() {
	if w == nil {
		return
	}
	w.mu.Lock()
	cancel := w.cancel
	w.cancel = nil
	w.mu.Unlock()
	if cancel != nil {
		cancel()
	}
}

// RunOnce 读取并补投暂停区间的记录，追平后交接；返回 nil 表示已恢复投递或无需继续。
func (w *Worker) RunOnce(ctx context.Context) error {
	if w == nil || w.intake == nil || w.store == nil {
		return nil
	}
	for {
		if err := ctx.Err(); err != nil {
			return err
		}
		if !w.needed() {
			return nil
		}
		paused, start := w.intake.PauseState()
		if !paused {
			return nil
		}
		st, err := w.stateFor(start)
		if err != nil {
			return err
		}

		flushedSeq, flushedOff := w.store.Flushed()
		if st.cursor.Seq < flushedSeq {
			before := st.cursor.Seq
			local := st.cursor
			records, after, scanErr := w.scan(&local, flushedOff, w.batchSize)
			if scanErr != nil {
				return scanErr
			}
			if len(records) == 0 {
				st.cursor = local
				if st.cursor.Seq == before {
					return errPending
				}
				continue
			}
			sent, replayErr := w.intake.Replay(ctx, records)
			if sent > 0 {
				st.cursor = after[sent-1]
			}
			if replayErr != nil {
				if errors.Is(replayErr, buffer.ErrNotPaused) {
					return nil
				}
				return replayErr
			}
			st.cursor = local
			continue
		}

		var next Cursor
		err = w.intake.Handover(ctx, w.store.Flush, func() ([]*recordv1.Record, uint64, error) {
			seq, off := w.store.Flushed()
			next = st.cursor
			records, _, scanErr := w.scan(&next, off, w.windowSize())
			if scanErr != nil {
				return nil, 0, scanErr
			}
			if next.Seq != seq {
				return nil, next.Seq, nil
			}
			return records, next.Seq, nil
		})
		switch {
		case err == nil:
			st.cursor = next
			w.state = state{}
			console.Debugf("backfill: resumed uploads after store seq %d", next.Seq)
			return nil
		case errors.Is(err, buffer.ErrNotPaused):
			return nil
		case errors.Is(err, buffer.ErrBackfillStale):
			// 交接窗口外仍有新提交，继续增量读取。
			continue
		case errors.Is(err, buffer.ErrBackfillBlocked):
			return errPending
		default:
			return err
		}
	}
}

// Complete 在收尾路径同步补完暂停区间；store 已关闭时不再刷盘。
func (w *Worker) Complete(ctx context.Context) error {
	if w == nil || w.intake == nil || w.store == nil {
		return nil
	}
	// 收尾路径独占游标；后台补投必须先停止。
	w.Stop()
	paused, start := w.intake.PauseState()
	if !paused {
		return nil
	}
	if w.degraded != nil && w.degraded.Load() {
		return fmt.Errorf("local store is degraded, records from storeSeq %d are not delivered", start)
	}
	if w.isSuspended != nil && w.isSuspended() {
		return fmt.Errorf("site is degraded, records from storeSeq %d are not delivered", start)
	}
	st, err := w.stateFor(start)
	if err != nil {
		return err
	}
	for {
		flushedSeq, flushedOff := w.store.Flushed()
		if st.cursor.Seq < flushedSeq {
			before := st.cursor.Seq
			local := st.cursor
			records, after, scanErr := w.scan(&local, flushedOff, w.batchSize)
			if scanErr != nil {
				return scanErr
			}
			if len(records) == 0 {
				st.cursor = local
				if st.cursor.Seq == before {
					return fmt.Errorf("local store ended before storeSeq %d", start)
				}
				continue
			}
			sent, replayErr := w.intake.Replay(ctx, records)
			if sent > 0 {
				st.cursor = after[sent-1]
			}
			if replayErr != nil {
				if errors.Is(replayErr, buffer.ErrNotPaused) {
					return nil
				}
				return fmt.Errorf("replay paused records: %w", replayErr)
			}
			st.cursor = local
			continue
		}

		var next Cursor
		err = w.intake.Handover(ctx, nil, func() ([]*recordv1.Record, uint64, error) {
			seq, off := w.store.Flushed()
			next = st.cursor
			records, _, scanErr := w.scan(&next, off, w.windowSize())
			if scanErr != nil {
				return nil, 0, scanErr
			}
			if next.Seq != seq {
				return nil, next.Seq, nil
			}
			return records, next.Seq, nil
		})
		switch {
		case err == nil:
			w.state = state{}
			console.Debugf("backfill: finished at store seq %d", next.Seq)
			return nil
		case errors.Is(err, buffer.ErrNotPaused):
			return nil
		case errors.Is(err, buffer.ErrBackfillStale):
			return errors.New("records were committed after the local store closed")
		case errors.Is(err, buffer.ErrBackfillBlocked):
			timer := time.NewTimer(pollInterval)
			select {
			case <-ctx.Done():
				timer.Stop()
				return ctx.Err()
			case <-timer.C:
			}
		default:
			return err
		}
	}
}

func (w *Worker) loop() {
	defer close(w.done)
	var retry time.Duration
	var fails int
	for {
		var timerC <-chan time.Time
		var timer *time.Timer
		if retry > 0 {
			timer = time.NewTimer(retry)
			timerC = timer.C
		}
		select {
		case <-w.ctx.Done():
			if timer != nil {
				timer.Stop()
			}
			return
		case <-w.stop:
			if timer != nil {
				timer.Stop()
			}
			return
		case <-w.wake:
		case <-timerC:
		}
		if timer != nil {
			timer.Stop()
		}

		ctx, ok := w.beginAttempt()
		if !ok {
			retry = 0
			continue
		}
		err := w.RunOnce(ctx)
		w.CancelAttempt()
		switch {
		case err == nil, errors.Is(err, context.Canceled), errors.Is(err, buffer.ErrNotPaused):
			retry, fails = 0, 0
		case errors.Is(err, errPending):
			retry = pollInterval
		case errors.Is(err, errCorrupt):
			// 确定性失败：停止重试并告警，避免只表现为长期暂停。
			console.Warningf("backfill: stopped after %v", err)
			retry, fails = 0, 0
		default:
			fails++
			retry = backoff(fails)
			console.Warningf("backfill: %v", err)
		}
	}
}

func (w *Worker) beginAttempt() (context.Context, bool) {
	w.mu.Lock()
	defer w.mu.Unlock()
	select {
	case <-w.stop:
		return nil, false
	default:
	}
	if !w.needed() {
		return nil, false
	}
	ctx, cancel := context.WithCancel(w.ctx)
	w.cancel = cancel
	return ctx, true
}

func (w *Worker) needed() bool {
	if w.intake == nil || w.store == nil {
		return false
	}
	if w.degraded != nil && w.degraded.Load() {
		return false
	}
	if w.isRunning != nil && !w.isRunning() {
		return false
	}
	if w.isSuspended != nil && w.isSuspended() {
		return false
	}
	paused, _ := w.intake.PauseState()
	return paused
}

func (w *Worker) stateFor(pauseSeq uint64) (*state, error) {
	if w.state.pauseSeq == pauseSeq {
		return &w.state, nil
	}
	cursor, err := w.newCursor(pauseSeq)
	if err != nil {
		return nil, err
	}
	w.state = state{pauseSeq: pauseSeq, cursor: cursor}
	console.Debugf("backfill: replaying from store seq %d", pauseSeq)
	return &w.state, nil
}

func (w *Worker) newCursor(target uint64) (Cursor, error) {
	if target == 0 {
		return Cursor{}, errors.New("backfill: storeSeq 0 is invalid")
	}
	flushedSeq, flushedOff := w.store.Flushed()
	if flushedSeq < target {
		return Cursor{}, errPending
	}
	cursor := Cursor{Offset: store.HeaderLen}
	if anchor, ok := w.store.Index().LocateFlushed(target-1, flushedOff); ok {
		// 锚点指向该序号记录本身的起点，从它前一条开始读。
		cursor = Cursor{Seq: anchor.StoreSeq - 1, Offset: anchor.Offset}
	}
	if skip := int64(target - 1 - cursor.Seq); skip > 0 {
		if _, _, err := w.scan(&cursor, flushedOff, int(skip)); err != nil {
			return Cursor{}, err
		}
		if cursor.Seq != target-1 {
			return Cursor{}, fmt.Errorf("backfill: store ended before storeSeq %d", target)
		}
	}
	return cursor, nil
}

// scan 从 cursor 顺序读取至多 want 条记录，跳过不可投递类型。
// after[i] 是产出第 i 条可投递记录后的游标，供部分投递后精确续读。
// 不会越过 limit 读取未刷盘的半条记录；cursor 只在前一条完整记录之后推进。
func (w *Worker) scan(
	cursor *Cursor,
	limit int64,
	want int,
) ([]*recordv1.Record, []Cursor, error) {
	if want <= 0 || cursor.Offset >= limit {
		return nil, nil, nil
	}
	reader, err := store.OpenOffset(w.runFile, cursor.Offset, limit)
	if err != nil {
		return nil, nil, err
	}
	defer func() { _ = reader.Close() }()

	var out []*recordv1.Record
	var after []Cursor
	for range want {
		payload, err := reader.Scan()
		switch {
		case errors.Is(err, io.EOF), errors.Is(err, store.ErrTruncated):
			return out, after, nil
		case errors.Is(err, store.ErrCorrupt):
			return out, after, fmt.Errorf("%w: %v", errCorrupt, err)
		case err != nil:
			return out, after, err
		}
		rec := &recordv1.Record{}
		if err := proto.Unmarshal(payload, rec); err != nil {
			return out, after, fmt.Errorf("%w: unmarshal record: %v", errCorrupt, err)
		}
		cursor.Seq++
		cursor.Offset = reader.Position()
		if w.onStored != nil {
			w.onStored(rec, cursor.Seq)
		}
		if record.IsImmediateUploadable(rec) {
			out = append(out, rec)
			after = append(after, *cursor)
		}
	}
	return out, after, nil
}

func (w *Worker) windowSize() int {
	return min(w.batchSize, windowMax)
}

func backoff(fails int) time.Duration {
	delay := retryBase
	for i := 1; i < fails && delay < retryMax; i++ {
		delay *= 2
	}
	return min(delay, retryMax)
}
