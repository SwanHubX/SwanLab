// 本文件承载预算暂停后的本地补投：从 pauseSeq 读取已落盘记录，按序投入上传通道，
// 追平后在 intake 锁内恢复实时投递。
package manager

import (
	"context"
	"errors"
	"fmt"
	"io"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

const (
	// backfillBatchMax 限制单次读取与补投的记录数，避免大 record_batch 占用过多内存。
	backfillBatchMax = 1024
	// backfillWindowMax 限制交接时持锁读取的尾部记录数。
	backfillWindowMax = 256
	// backfillPoll 是等待刷盘推进或上传槽位释放的轮询间隔。
	backfillPoll = 50 * time.Millisecond
	// backfillRetry / backfillRetryMax 是失败重试的初始与最大退避。
	backfillRetry    = 200 * time.Millisecond
	backfillRetryMax = 2 * time.Second
)

var (
	// errBackfillPending 表示本次没有可推进的工作，稍后重试即可。
	errBackfillPending = errors.New("manager: backfill pending")
	// errBackfillCorrupt 表示本地记录无法解析，重试不会成功。
	errBackfillCorrupt = errors.New("manager: backfill stopped on corrupt store")
)

// backfillCursor 指向下一条待读记录：seq 为已读取的最后一条 storeSeq。
type backfillCursor struct {
	seq    uint64
	offset int64
}

// backfillState 只在补投协程或已停止补投的收尾路径中访问。
type backfillState struct {
	pauseSeq uint64
	cursor   backfillCursor
}

// startBackfill 启动补投协程；只在在线会话装配完成后调用一次。
func (s *Session) startBackfill() {
	s.backfillWake = make(chan struct{}, 1)
	s.backfillStop = make(chan struct{})
	s.backfillDone = make(chan struct{})
	go s.backfillLoop()
}

// stopBackfill 取消当前尝试并等待协程退出；可重复调用。
func (s *Session) stopBackfill() {
	if s.backfillStop == nil {
		return
	}
	s.backfillOnce.Do(func() { close(s.backfillStop) })
	s.cancelBackfillAttempt()
	<-s.backfillDone
}

// backfillLoop 在暂停、健康与刷盘条件变化时重试补投。
func (s *Session) backfillLoop() {
	defer close(s.backfillDone)
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
		case <-s.ctx.Done():
			if timer != nil {
				timer.Stop()
			}
			return
		case <-s.backfillStop:
			if timer != nil {
				timer.Stop()
			}
			return
		case <-s.backfillWake:
		case <-timerC:
		}
		if timer != nil {
			timer.Stop()
		}

		ctx, ok := s.beginBackfillAttempt()
		if !ok {
			retry = 0
			continue
		}
		err := s.runBackfill(ctx)
		s.cancelBackfillAttempt()
		switch {
		case err == nil, errors.Is(err, context.Canceled), errors.Is(err, buffer.ErrNotPaused):
			retry, fails = 0, 0
		case errors.Is(err, errBackfillPending):
			retry = backfillPoll
		case errors.Is(err, errBackfillCorrupt):
			// 确定性失败：停止重试并告警，避免只表现为长期暂停。
			console.Warningf("backfill: stopped after %v", err)
			retry, fails = 0, 0
		default:
			fails++
			retry = backfillBackoff(fails)
			console.Warningf("backfill: %v", err)
		}
	}
}

// beginBackfillAttempt 在确认补投仍有必要时创建可取消的尝试 ctx。
func (s *Session) beginBackfillAttempt() (context.Context, bool) {
	s.backfillMu.Lock()
	defer s.backfillMu.Unlock()
	select {
	case <-s.backfillStop:
		return nil, false
	default:
	}
	if !s.backfillNeeded() {
		return nil, false
	}
	ctx, cancel := context.WithCancel(s.ctx)
	s.backfillCancel = cancel
	return ctx, true
}

// cancelBackfillAttempt 取消并清除当前尝试；补投结束与健康变化时调用。
func (s *Session) cancelBackfillAttempt() {
	s.backfillMu.Lock()
	cancel := s.backfillCancel
	s.backfillCancel = nil
	s.backfillMu.Unlock()
	if cancel != nil {
		cancel()
	}
}

// wakeBackfill 唤醒补投协程；可在任意协程调用。
func (s *Session) wakeBackfill() {
	if s.backfillWake == nil {
		return
	}
	select {
	case s.backfillWake <- struct{}{}:
	default:
	}
}

// backfillNeeded 判断当前是否应继续补投。
func (s *Session) backfillNeeded() bool {
	if s.intake == nil || s.store == nil || s.uploader == nil {
		return false
	}
	if s.cfg.SkipStore || s.storeDegraded.Load() {
		return false
	}
	if s.stateSnapshot() != RunStateRunning || s.shouldSuspendUploads() {
		return false
	}
	paused, _ := s.intake.PauseState()
	return paused
}

// backfillBackoff 返回连续失败后的重试间隔。
func backfillBackoff(fails int) time.Duration {
	delay := backfillRetry
	for i := 1; i < fails && delay < backfillRetryMax; i++ {
		delay *= 2
	}
	return min(delay, backfillRetryMax)
}

// runBackfill 读取并补投暂停区间的记录，追平后交接；返回 nil 表示已恢复投递。
func (s *Session) runBackfill(ctx context.Context) error {
	for {
		if err := ctx.Err(); err != nil {
			return err
		}
		if !s.backfillNeeded() {
			return nil
		}
		paused, start := s.intake.PauseState()
		if !paused {
			return nil
		}
		state, err := s.backfillStateFor(start)
		if err != nil {
			return err
		}

		flushedSeq, flushedOff := s.store.Flushed()
		if state.cursor.seq < flushedSeq {
			before := state.cursor.seq
			local := state.cursor
			records, after, scanErr := s.scanBackfill(&local, flushedOff, s.backfillBatchSize())
			if scanErr != nil {
				return scanErr
			}
			if len(records) == 0 {
				state.cursor = local
				if state.cursor.seq == before {
					return errBackfillPending
				}
				continue
			}
			sent, replayErr := s.intake.Replay(ctx, records)
			if sent > 0 {
				state.cursor = after[sent-1]
			}
			if replayErr != nil {
				if errors.Is(replayErr, buffer.ErrNotPaused) {
					return nil
				}
				return replayErr
			}
			state.cursor = local
			continue
		}

		var next backfillCursor
		err = s.intake.Handover(ctx, s.store.Flush, func() ([]*recordv1.Record, uint64, error) {
			seq, off := s.store.Flushed()
			next = state.cursor
			records, _, scanErr := s.scanBackfill(&next, off, s.backfillWindow())
			if scanErr != nil {
				return nil, 0, scanErr
			}
			if next.seq != seq {
				return nil, next.seq, nil
			}
			return records, next.seq, nil
		})
		switch {
		case err == nil:
			state.cursor = next
			s.backfillState = backfillState{}
			console.Debugf("backfill: resumed uploads after store seq %d", next.seq)
			return nil
		case errors.Is(err, buffer.ErrNotPaused):
			return nil
		case errors.Is(err, buffer.ErrBackfillStale):
			// 交接窗口外仍有新提交，继续增量读取。
			continue
		case errors.Is(err, buffer.ErrBackfillBlocked):
			return errBackfillPending
		default:
			return err
		}
	}
}

// completeBackfill 在收尾路径同步补完暂停区间；store 已关闭时不再刷盘。
func (s *Session) completeBackfill(ctx context.Context) error {
	if s.intake == nil || s.store == nil || s.uploader == nil {
		return nil
	}
	// 收尾路径独占游标；后台补投必须先停止。
	s.stopBackfill()
	paused, start := s.intake.PauseState()
	if !paused {
		return nil
	}
	if s.storeDegraded.Load() {
		return fmt.Errorf("%w: local store is degraded, records from storeSeq %d are not delivered", ErrRunFailed, start)
	}
	if s.shouldSuspendUploads() {
		return fmt.Errorf("%w: site is degraded, records from storeSeq %d are not delivered", ErrRunFailed, start)
	}
	state, err := s.backfillStateFor(start)
	if err != nil {
		return fmt.Errorf("%w: %v", ErrRunFailed, err)
	}
	for {
		flushedSeq, flushedOff := s.store.Flushed()
		if state.cursor.seq < flushedSeq {
			before := state.cursor.seq
			local := state.cursor
			records, after, scanErr := s.scanBackfill(&local, flushedOff, s.backfillBatchSize())
			if scanErr != nil {
				return fmt.Errorf("%w: %v", ErrRunFailed, scanErr)
			}
			if len(records) == 0 {
				state.cursor = local
				if state.cursor.seq == before {
					return fmt.Errorf("%w: local store ended before storeSeq %d", ErrRunFailed, start)
				}
				continue
			}
			sent, replayErr := s.intake.Replay(ctx, records)
			if sent > 0 {
				state.cursor = after[sent-1]
			}
			if replayErr != nil {
				if errors.Is(replayErr, buffer.ErrNotPaused) {
					return nil
				}
				return fmt.Errorf("%w: replay paused records: %v", ErrRunFailed, replayErr)
			}
			state.cursor = local
			continue
		}

		var next backfillCursor
		err = s.intake.Handover(ctx, nil, func() ([]*recordv1.Record, uint64, error) {
			seq, off := s.store.Flushed()
			next = state.cursor
			records, _, scanErr := s.scanBackfill(&next, off, s.backfillWindow())
			if scanErr != nil {
				return nil, 0, scanErr
			}
			if next.seq != seq {
				return nil, next.seq, nil
			}
			return records, next.seq, nil
		})
		switch {
		case err == nil:
			s.backfillState = backfillState{}
			console.Debugf("backfill: finished at store seq %d", next.seq)
			return nil
		case errors.Is(err, buffer.ErrNotPaused):
			return nil
		case errors.Is(err, buffer.ErrBackfillStale):
			return fmt.Errorf("%w: records were committed after the local store closed", ErrRunFailed)
		case errors.Is(err, buffer.ErrBackfillBlocked):
			timer := time.NewTimer(backfillPoll)
			select {
			case <-ctx.Done():
				timer.Stop()
				return ctx.Err()
			case <-timer.C:
			}
		default:
			return fmt.Errorf("%w: %v", ErrRunFailed, err)
		}
	}
}

// backfillStateFor 返回当前暂停区间的游标；区间变化时重新定位。
func (s *Session) backfillStateFor(pauseSeq uint64) (*backfillState, error) {
	if s.backfillState.pauseSeq == pauseSeq {
		return &s.backfillState, nil
	}
	cursor, err := s.newBackfillCursor(pauseSeq)
	if err != nil {
		return nil, err
	}
	s.backfillState = backfillState{pauseSeq: pauseSeq, cursor: cursor}
	console.Debugf("backfill: replaying from store seq %d", pauseSeq)
	return &s.backfillState, nil
}

// newBackfillCursor 定位到 storeSeq 对应记录的起点；记录尚未刷盘时返回 pending。
func (s *Session) newBackfillCursor(target uint64) (backfillCursor, error) {
	if target == 0 {
		return backfillCursor{}, errors.New("backfill: storeSeq 0 is invalid")
	}
	flushedSeq, flushedOff := s.store.Flushed()
	if flushedSeq < target {
		return backfillCursor{}, errBackfillPending
	}
	cursor := backfillCursor{offset: store.HeaderLen}
	if anchor, ok := s.store.Index().LocateFlushed(target-1, flushedOff); ok {
		// 锚点指向该序号记录本身的起点，从它前一条开始读。
		cursor = backfillCursor{seq: anchor.StoreSeq - 1, offset: anchor.Offset}
	}
	if skip := int64(target - 1 - cursor.seq); skip > 0 {
		if _, _, err := s.scanBackfill(&cursor, flushedOff, int(skip)); err != nil {
			return backfillCursor{}, err
		}
		if cursor.seq != target-1 {
			return backfillCursor{}, fmt.Errorf("backfill: store ended before storeSeq %d", target)
		}
	}
	return cursor, nil
}

// scanBackfill 从 cursor 顺序读取至多 max 条记录，跳过不可投递类型。
// after[i] 是产出第 i 条可投递记录后的游标，供部分投递后精确续读。
// 不会越过 limit 读取未刷盘的半条记录；cursor 只在前一条完整记录之后推进。
func (s *Session) scanBackfill(
	cursor *backfillCursor,
	limit int64,
	want int,
) ([]*recordv1.Record, []backfillCursor, error) {
	if want <= 0 || cursor.offset >= limit {
		return nil, nil, nil
	}
	reader, err := store.OpenOffset(s.runFile, cursor.offset, limit)
	if err != nil {
		return nil, nil, err
	}
	defer func() { _ = reader.Close() }()

	var out []*recordv1.Record
	var after []backfillCursor
	for range want {
		payload, err := reader.Scan()
		switch {
		case errors.Is(err, io.EOF), errors.Is(err, store.ErrTruncated):
			return out, after, nil
		case errors.Is(err, store.ErrCorrupt):
			return out, after, fmt.Errorf("%w: %v", errBackfillCorrupt, err)
		case err != nil:
			return out, after, err
		}
		rec := &recordv1.Record{}
		if err := proto.Unmarshal(payload, rec); err != nil {
			return out, after, fmt.Errorf("%w: unmarshal record: %v", errBackfillCorrupt, err)
		}
		cursor.seq++
		cursor.offset = reader.Position()
		if backfillUploadable(rec) {
			out = append(out, rec)
			after = append(after, *cursor)
		}
	}
	return out, after, nil
}

// backfillUploadable 过滤只落本地与延迟投递的记录。
func backfillUploadable(rec *recordv1.Record) bool {
	switch record.TypeOf(rec) {
	case record.TypeStart, record.TypeFinish:
		return false
	case record.TypeSave:
		// END 保存在 Close 时由 intake 统一投递，补投跳过以免重复。
		return !deferredSave(rec)
	default:
		return true
	}
}

// backfillBatchSize 返回单次补投的记录数。
func (s *Session) backfillBatchSize() int {
	batch := int(s.cfg.RecordBatch)
	if batch <= 0 {
		batch = backfillBatchMax
	}
	return min(batch, backfillBatchMax)
}

// backfillWindow 返回交接时持锁读取的最大记录数。
func (s *Session) backfillWindow() int {
	return min(s.backfillBatchSize(), backfillWindowMax)
}
