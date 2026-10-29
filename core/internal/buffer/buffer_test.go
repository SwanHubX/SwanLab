package buffer

import (
	"context"
	"errors"
	"fmt"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func scalarRecord(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Scalar{Scalar: &datav1.ScalarRecord{
			Key:   "train/loss",
			Step:  num,
			Type:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
			Value: &datav1.ScalarValue{Number: float64(num)},
		}},
	}
}

func startRecord() *recordv1.Record {
	return &recordv1.Record{Num: record.NumStart, RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}}}
}

func finishRecord() *recordv1.Record {
	return &recordv1.Record{Num: record.NumFinish, RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{}}}
}

func logRecord(num int64) *recordv1.Record {
	return &recordv1.Record{Num: num, RecordType: &recordv1.Record_Log{Log: &terminalv1.LogRecord{Line: "line", Epoch: num}}}
}

func saveRecord(num int64, policy savev1.SavePolicy) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{Save: &savev1.SaveRecord{
			Name:   fmt.Sprintf("file-%d", num),
			Policy: policy,
			Type:   savev1.SaveType_SAVE_TYPE_CUSTOM,
		}},
	}
}

func demandFor(recs []*recordv1.Record) Demand {
	var bytes int64
	for _, rec := range recs {
		bytes += record.Size(rec)
	}
	return Demand{Store: len(recs), Upload: len(recs), Bytes: bytes}
}

type storeSink struct {
	ch     chan store.Request
	mu     sync.Mutex
	got    []store.Request
	notify chan struct{}
	done   chan struct{}
}

func newStoreSink(size int) *storeSink {
	s := &storeSink{ch: make(chan store.Request, size), notify: make(chan struct{}, 1), done: make(chan struct{})}
	go func() {
		defer close(s.done)
		for req := range s.ch {
			s.mu.Lock()
			s.got = append(s.got, req)
			s.mu.Unlock()
			select {
			case s.notify <- struct{}{}:
			default:
			}
		}
	}()
	return s
}

func (s *storeSink) stop() {
	close(s.ch)
	<-s.done
}

func (s *storeSink) snapshot() []store.Request {
	s.mu.Lock()
	defer s.mu.Unlock()
	return append([]store.Request(nil), s.got...)
}

func newIntake(t *testing.T, opts Options) *Intake {
	t.Helper()
	in, err := New(opts)
	if err != nil {
		t.Fatalf("New: %v", err)
	}
	return in
}

func drainUpload(in *Intake) []*recordv1.Record {
	var out []*recordv1.Record
	for {
		select {
		case rec, ok := <-in.UploadOut():
			if !ok {
				return out
			}
			out = append(out, rec)
		default:
			return out
		}
	}
}

func waitStore(t *testing.T, sink *storeSink, want int) []store.Request {
	t.Helper()
	deadline := time.After(time.Second)
	for len(sink.snapshot()) < want {
		select {
		case <-sink.notify:
		case <-deadline:
			t.Fatalf("store received %d records, want %d", len(sink.snapshot()), want)
		}
	}
	return sink.snapshot()
}

func TestNewValidatesOptions(t *testing.T) {
	if _, err := New(Options{}); err == nil {
		t.Fatal("New without budget returned nil error")
	}
	b := budget.New(budget.Limits{})
	if _, err := New(Options{Budget: b}); err == nil {
		t.Fatal("New without store channel returned nil error")
	}
	if _, err := New(Options{Budget: b, SkipStore: true}); err != nil {
		t.Fatalf("New with skip_store: %v", err)
	}
}

func TestCommitEmptyAndNil(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 4})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	tk, err := in.Reserve(context.Background(), Demand{Store: 1, Upload: 1, Bytes: 10})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err = in.Commit(tk, nil); err != nil {
		t.Fatalf("Commit empty: %v", err)
	}
	if records, bytes := b.Usage(); records != 0 || bytes != 0 {
		t.Fatalf("budget usage = (%d, %d), want (0, 0)", records, bytes)
	}

	tk, err = in.Reserve(context.Background(), Demand{Store: 1, Upload: 1, Bytes: 10})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err = in.Commit(tk, []*recordv1.Record{nil}); err == nil {
		t.Fatal("Commit with nil record returned nil error")
	}
	if records, bytes := b.Usage(); records != 0 || bytes != 0 {
		t.Fatalf("budget usage after nil record = (%d, %d), want (0, 0)", records, bytes)
	}
}

func TestCommitAfterClose(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 4})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	tk, err := in.Reserve(context.Background(), Demand{Store: 1, Upload: 1, Bytes: 10})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if err := in.Commit(tk, []*recordv1.Record{scalarRecord(1)}); !errors.Is(err, ErrClosed) {
		t.Fatalf("Commit err = %v, want ErrClosed", err)
	}
	if records, bytes := b.Usage(); records != 0 || bytes != 0 {
		t.Fatalf("budget usage = (%d, %d), want (0, 0)", records, bytes)
	}
}

func TestCloseRetriesPendingEndSaves(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 1})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	first := []*recordv1.Record{scalarRecord(1)}
	tk, err := in.Reserve(context.Background(), demandFor(first))
	if err != nil {
		t.Fatalf("Reserve first: %v", err)
	}
	if err = in.Commit(tk, first); err != nil {
		t.Fatalf("Commit first: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload got %d, want 1", got)
	}

	// 预算已满 → 暂停投递；END 保存落本地并暂存为未扣配额。
	save := []*recordv1.Record{saveRecord(2, savev1.SavePolicy_SAVE_POLICY_END)}
	tk, err = in.Reserve(context.Background(), demandFor(save))
	if err != nil {
		t.Fatalf("Reserve save: %v", err)
	}
	if err = in.Commit(tk, save); err != nil {
		t.Fatalf("Commit save: %v", err)
	}
	if got := waitStore(t, sink, 2); len(got) != 2 {
		t.Fatalf("store got %d, want 2", len(got))
	}

	// 预算仍被第一条记录占用：Close 有界等待超时，暂存记录放回队列。
	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
	defer cancel()
	if err := in.Close(ctx); !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("Close err = %v, want context.DeadlineExceeded", err)
	}

	// 模拟 uploader 确认第一条后归还配额，重试 Close 应补扣并投递 END 保存。
	records, bytes := b.Usage()
	b.Release(records, bytes)
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("second Close: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload after close got %d, want 1", got)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("budget usage = %d, want 1 (deferred save charged at close)", records)
	}
}

func TestCommitUploadOnlyWaitsForBudget(t *testing.T) {
	b := budget.New(budget.Limits{Records: 1})
	if !b.TryAcquire(1, 0) {
		t.Fatal("seed TryAcquire = false")
	}
	in := newIntake(t, Options{SkipStore: true, Budget: b, UploadQueue: 2})

	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
	defer cancel()
	if err := in.CommitUploadOnly(ctx, logRecord(1)); !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("CommitUploadOnly err = %v, want context.DeadlineExceeded", err)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("budget usage = %d, want 1", records)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload got %d, want 0", got)
	}
}

func TestCommitChargesActualBytes(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{Bytes: 1 << 20})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	rec := scalarRecord(1)
	// 故意低估预留字节，Commit 必须按实际序列化大小补扣。
	tk, err := in.Reserve(context.Background(), Demand{Store: 1, Upload: 1, Bytes: 1})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := in.Commit(tk, []*recordv1.Record{rec}); err != nil {
		t.Fatalf("Commit: %v", err)
	}
	_, bytes := b.Usage()
	if want := record.Size(rec); bytes != want {
		t.Fatalf("budget bytes = %d, want %d", bytes, want)
	}
}

func TestPendingEndSaveUploadsAtClose(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	save := []*recordv1.Record{saveRecord(1, savev1.SavePolicy_SAVE_POLICY_END)}
	tk, err := in.Reserve(context.Background(), demandFor(save))
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := in.Commit(tk, save); err != nil {
		t.Fatalf("Commit: %v", err)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload before close got %d, want 0", got)
	}
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload after close got %d, want 1", got)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("budget usage = %d, want 1", records)
	}
}

func TestCommitRoutes(t *testing.T) {
	tests := []struct {
		name       string
		records    []*recordv1.Record
		demand     Demand
		wantStore  int
		wantUpload int
		flush      bool // Close 后再统计上传通道
	}{
		{
			name:       "scalar both channels",
			records:    []*recordv1.Record{scalarRecord(1)},
			wantStore:  1,
			wantUpload: 1,
		},
		{
			name:      "start local only",
			records:   []*recordv1.Record{startRecord()},
			demand:    Demand{Store: 1},
			wantStore: 1,
		},
		{
			name:      "finish local only",
			records:   []*recordv1.Record{finishRecord()},
			demand:    Demand{Store: 1},
			wantStore: 1,
		},
		{
			name:       "log both channels",
			records:    []*recordv1.Record{logRecord(5)},
			wantStore:  1,
			wantUpload: 1,
		},
		{
			name:       "save now both channels",
			records:    []*recordv1.Record{saveRecord(6, savev1.SavePolicy_SAVE_POLICY_NOW)},
			wantStore:  1,
			wantUpload: 1,
		},
		{
			name:      "save end deferred until close",
			records:   []*recordv1.Record{saveRecord(7, savev1.SavePolicy_SAVE_POLICY_END)},
			wantStore: 1,
			flush:     true,
			// 延迟上传仍占用上传预算，Close 后才进入通道。
			wantUpload: 1,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			sink := newStoreSink(8)
			defer sink.stop()
			b := budget.New(budget.Limits{})
			in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 8})

			d := tc.demand
			if len(tc.records) > 0 && (d.Store == 0 && d.Upload == 0) {
				d = demandFor(tc.records)
			}
			tk, err := in.Reserve(context.Background(), d)
			if err != nil {
				t.Fatalf("Reserve: %v", err)
			}
			if err := in.Commit(tk, tc.records); err != nil {
				t.Fatalf("Commit: %v", err)
			}

			if got := waitStore(t, sink, tc.wantStore); len(got) != tc.wantStore {
				t.Fatalf("store got %d, want %d", len(got), tc.wantStore)
			}
			if !tc.flush {
				if got := len(drainUpload(in)); got != tc.wantUpload {
					t.Fatalf("upload got %d, want %d", got, tc.wantUpload)
				}
				return
			}
			if got := len(drainUpload(in)); got != 0 {
				t.Fatalf("upload before close got %d, want 0", got)
			}
			if err := in.Close(context.Background()); err != nil {
				t.Fatalf("Close: %v", err)
			}
			if got := len(drainUpload(in)); got != tc.wantUpload {
				t.Fatalf("upload after close got %d, want %d", got, tc.wantUpload)
			}
		})
	}
}

func TestSkipStoreRoutesDataToUpload(t *testing.T) {
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{SkipStore: true, Budget: b, UploadQueue: 4})
	recs := []*recordv1.Record{scalarRecord(1), scalarRecord(2)}
	tk, err := in.Reserve(context.Background(), Demand{Upload: 2, Bytes: 32})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err = in.Commit(tk, recs); err != nil {
		t.Fatalf("Commit: %v", err)
	}
	if got := len(drainUpload(in)); got != 2 {
		t.Fatalf("upload got %d, want 2", got)
	}
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

func TestPauseAndResume(t *testing.T) {
	sink := newStoreSink(16)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 2})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 8})

	first := []*recordv1.Record{scalarRecord(1), scalarRecord(2)}
	tk, err := in.Reserve(context.Background(), demandFor(first))
	if err != nil {
		t.Fatalf("Reserve first: %v", err)
	}
	if err = in.Commit(tk, first); err != nil {
		t.Fatalf("Commit first: %v", err)
	}

	second := []*recordv1.Record{scalarRecord(3), scalarRecord(4)}
	start := time.Now()
	tk, err = in.Reserve(context.Background(), demandFor(second))
	if err != nil {
		t.Fatalf("Reserve second: %v", err)
	}
	if elapsed := time.Since(start); elapsed > time.Second {
		t.Fatalf("Reserve blocked for %v while a local copy exists", elapsed)
	}
	if err = in.Commit(tk, second); err != nil {
		t.Fatalf("Commit second: %v", err)
	}

	paused, seq := in.PauseState()
	if !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}
	if got := waitStore(t, sink, 4); len(got) != 4 {
		t.Fatalf("store got %d, want 4", len(got))
	}
	if got := len(drainUpload(in)); got != 2 {
		t.Fatalf("upload got %d, want 2 (paused batch must not upload)", got)
	}
	if records, _ := b.Usage(); records != 2 {
		t.Fatalf("budget usage = %d, want 2", records)
	}

	// 模拟上传确认后归还配额，再恢复投递。
	records, bytes := b.Usage()
	b.Release(records, bytes)
	in.ResumeUpload()
	if paused, seq := in.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState after resume = (%v, %d), want (false, 0)", paused, seq)
	}

	third := []*recordv1.Record{scalarRecord(5)}
	tk, err = in.Reserve(context.Background(), demandFor(third))
	if err != nil {
		t.Fatalf("Reserve third: %v", err)
	}
	if err := in.Commit(tk, third); err != nil {
		t.Fatalf("Commit third: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload after resume got %d, want 1", got)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("budget usage = %d, want 1", records)
	}
}

func TestCommitAfterResumeKeepsPausedTicket(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 1})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	first := []*recordv1.Record{scalarRecord(1)}
	tk, err := in.Reserve(context.Background(), demandFor(first))
	if err != nil {
		t.Fatalf("Reserve first: %v", err)
	}
	if err = in.Commit(tk, first); err != nil {
		t.Fatalf("Commit first: %v", err)
	}

	second := []*recordv1.Record{scalarRecord(2)}
	tk, err = in.Reserve(context.Background(), demandFor(second))
	if err != nil {
		t.Fatalf("Reserve second: %v", err)
	}
	// 模拟恢复投递发生在 Commit 之前：票据仍按暂停处理，不能报超额。
	in.ResumeUpload()
	if err = in.Commit(tk, second); err != nil {
		t.Fatalf("Commit after resume: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload got %d, want 1", got)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("budget usage = %d, want 1", records)
	}
	if got := waitStore(t, sink, 2); len(got) != 2 {
		t.Fatalf("store got %d, want 2", len(got))
	}
	if paused, seq := in.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState after resume = (%v, %d), want (false, 0)", paused, seq)
	}
}

func TestDegradedUsesBoundedBackpressure(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 1})
	if !b.TryAcquire(1, 0) {
		t.Fatal("seed TryAcquire = false")
	}
	var degraded atomic.Bool
	degraded.Store(true)
	in := newIntake(t, Options{
		Store:       sink.ch,
		Budget:      b,
		Degraded:    &degraded,
		UploadQueue: 4,
		QueueWait:   30 * time.Millisecond,
	})

	start := time.Now()
	_, err := in.Reserve(context.Background(), Demand{Store: 1, Upload: 1, Bytes: 10})
	elapsed := time.Since(start)
	if !errors.Is(err, ErrOverloaded) {
		t.Fatalf("Reserve err = %v, want ErrOverloaded", err)
	}
	if elapsed < 20*time.Millisecond {
		t.Fatalf("Reserve returned in %v, want bounded wait", elapsed)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("budget usage = %d, want 1 (failed reserve must not change state)", records)
	}
	if got := len(sink.snapshot()); got != 0 {
		t.Fatalf("store received %d records, want 0", got)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload got %d records, want 0", got)
	}
}

func TestStoreQueueFullBackpressuresWithoutLoss(t *testing.T) {
	storeCh := make(chan store.Request, 4)
	for i := 0; i < cap(storeCh); i++ {
		storeCh <- store.Request{StoreSeq: uint64(i + 1), Payload: []byte("x")}
	}
	var degraded atomic.Bool
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{
		Store:             storeCh,
		Budget:            b,
		Degraded:          &degraded,
		UploadQueue:       8,
		QueueWait:         20 * time.Millisecond,
		StorePollInterval: time.Millisecond,
	})

	recs := []*recordv1.Record{scalarRecord(1), scalarRecord(2)}
	start := time.Now()
	_, err := in.Reserve(context.Background(), demandFor(recs))
	if !errors.Is(err, ErrOverloaded) {
		t.Fatalf("Reserve err = %v, want ErrOverloaded", err)
	}
	if elapsed := time.Since(start); elapsed < 15*time.Millisecond {
		t.Fatalf("Reserve returned in %v, want bounded wait", elapsed)
	}
	if degraded.Load() {
		t.Fatal("degraded = true after capacity pressure, want false (no loss/降级)")
	}
	if records, bytes := b.Usage(); records != 0 || bytes != 0 {
		t.Fatalf("budget usage = (%d, %d), want (0, 0)", records, bytes)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload got %d records, want 0", got)
	}
	if got := len(storeCh); got != 4 {
		t.Fatalf("store queue length = %d, want 4 (no writes)", got)
	}
}

func TestBatchLargerThanStoreQueueIsRejected(t *testing.T) {
	storeCh := make(chan store.Request, 2)
	var degraded atomic.Bool
	in := newIntake(t, Options{
		Store:       storeCh,
		Budget:      budget.New(budget.Limits{}),
		Degraded:    &degraded,
		UploadQueue: 4,
	})
	_, err := in.Reserve(context.Background(), Demand{Store: 3, Upload: 3, Bytes: 30})
	if !errors.Is(err, ErrBatchTooLarge) {
		t.Fatalf("Reserve err = %v, want ErrBatchTooLarge", err)
	}
	if degraded.Load() {
		t.Fatal("degraded = true, want false")
	}
	if got := len(storeCh); got != 0 {
		t.Fatalf("store queue length = %d, want 0", got)
	}
}

func TestPausedBatchFallsBackToUploadWhenStoreDegrades(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 1})
	var degraded atomic.Bool
	in := newIntake(t, Options{Store: sink.ch, Budget: b, Degraded: &degraded, UploadQueue: 4})

	first := []*recordv1.Record{scalarRecord(1)}
	tk, err := in.Reserve(context.Background(), demandFor(first))
	if err != nil {
		t.Fatalf("Reserve first: %v", err)
	}
	if err = in.Commit(tk, first); err != nil {
		t.Fatalf("Commit first: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload got %d, want 1", got)
	}

	// 预算满 → 第二批按暂停处理；随后本地存储降级，Commit 必须改走上传。
	second := []*recordv1.Record{scalarRecord(2)}
	tk, err = in.Reserve(context.Background(), demandFor(second))
	if err != nil {
		t.Fatalf("Reserve second: %v", err)
	}
	degraded.Store(true)
	if err = in.Commit(tk, second); err != nil {
		t.Fatalf("Commit second: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload after degrade got %d, want 1 (paused batch must not be dropped)", got)
	}
	if got := waitStore(t, sink, 1); len(got) != 1 {
		t.Fatalf("store got %d, want 1 (only the first record)", len(got))
	}
	if records, _ := b.Usage(); records != 2 {
		t.Fatalf("budget usage = %d, want 2 (temporary overrun is allowed to avoid loss)", records)
	}
}

func TestDegradedDoesNotMutateSkipStore(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{})
	var degraded atomic.Bool
	in := newIntake(t, Options{Store: sink.ch, Budget: b, Degraded: &degraded, UploadQueue: 4})
	if in.skipStore {
		t.Fatal("skipStore = true before degrade, want false")
	}

	// 模拟磁盘写满：只有运行期 storeDegraded 置位，配置保持原样。
	degraded.Store(true)
	recs := []*recordv1.Record{scalarRecord(1)}
	tk, err := in.Reserve(context.Background(), demandFor(recs))
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err = in.Commit(tk, recs); err != nil {
		t.Fatalf("Commit: %v", err)
	}
	if in.skipStore {
		t.Fatal("storeDegraded must not rewrite the skip_store setting")
	}
	if in.storeOut == nil {
		t.Fatal("store channel must stay configured after degrade")
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload got %d, want 1", got)
	}
	if got := len(sink.snapshot()); got != 0 {
		t.Fatalf("store received %d records, want 0", got)
	}
}

func TestCommitUploadOnly(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	rec := logRecord(9)
	if err := in.CommitUploadOnly(context.Background(), rec); err != nil {
		t.Fatalf("CommitUploadOnly: %v", err)
	}
	if got := len(drainUpload(in)); got != 1 {
		t.Fatalf("upload got %d, want 1", got)
	}
	if got := len(sink.snapshot()); got != 0 {
		t.Fatalf("store received %d records, want 0", got)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("budget usage = %d, want 1", records)
	}
}

func TestCommitUploadOnlyAfterClose(t *testing.T) {
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{SkipStore: true, Budget: b, UploadQueue: 2})
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}
	err := in.CommitUploadOnly(context.Background(), logRecord(1))
	if !errors.Is(err, ErrClosed) {
		t.Fatalf("CommitUploadOnly err = %v, want ErrClosed", err)
	}
	if records, _ := b.Usage(); records != 0 {
		t.Fatalf("budget usage = %d, want 0 after rejected commit", records)
	}
}

func TestCloseIsIdempotentAndRejectsCommits(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("second Close: %v", err)
	}
	if _, err := in.Reserve(context.Background(), Demand{Store: 1, Upload: 1}); !errors.Is(err, ErrClosed) {
		t.Fatalf("Reserve after close err = %v, want ErrClosed", err)
	}
}

func TestOverReservedCommitHasNoSideEffects(t *testing.T) {
	sink := newStoreSink(4)
	defer sink.stop()
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	tk, err := in.Reserve(context.Background(), Demand{Store: 1, Upload: 1, Bytes: 10})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	err = in.Commit(tk, []*recordv1.Record{scalarRecord(1), scalarRecord(2)})
	if !errors.Is(err, ErrOverReserved) {
		t.Fatalf("Commit err = %v, want ErrOverReserved", err)
	}
	if got := len(sink.snapshot()); got != 0 {
		t.Fatalf("store received %d records, want 0", got)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload got %d records, want 0", got)
	}
	if records, bytes := b.Usage(); records != 0 || bytes != 0 {
		t.Fatalf("budget usage = (%d, %d), want (0, 0)", records, bytes)
	}
}

func TestReleaseUnusedReservation(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 4})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	recs := []*recordv1.Record{scalarRecord(1), scalarRecord(2)}
	tk, err := in.Reserve(context.Background(), Demand{Store: 3, Upload: 3, Bytes: 30})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err = in.Commit(tk, recs); err != nil {
		t.Fatalf("Commit: %v", err)
	}
	if records, _ := b.Usage(); records != 2 {
		t.Fatalf("budget usage = %d, want 2", records)
	}
	// 预留的 1 条本地空间也必须归还，否则后续 Reserve 会误判队列已满。
	tk2, err := in.Reserve(context.Background(), Demand{Store: 2, Upload: 2, Bytes: 20})
	if err != nil {
		t.Fatalf("second Reserve: %v", err)
	}
	in.Release(tk2)
	if got := waitStore(t, sink, 2); len(got) != 2 {
		t.Fatalf("store got %d, want 2", len(got))
	}
}

func TestConcurrentCommits(t *testing.T) {
	const workers = 8
	sink := newStoreSink(workers)
	defer sink.stop()
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: workers})

	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			recs := []*recordv1.Record{scalarRecord(int64(i + 1))}
			tk, err := in.Reserve(context.Background(), demandFor(recs))
			if err != nil {
				t.Errorf("Reserve: %v", err)
				return
			}
			if err = in.Commit(tk, recs); err != nil {
				t.Errorf("Commit: %v", err)
			}
		}(i)
	}
	wg.Wait()
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}

	got := waitStore(t, sink, workers)
	if len(got) != workers {
		t.Fatalf("store got %d, want %d", len(got), workers)
	}
	seen := make(map[uint64]bool, workers)
	for _, req := range got {
		if req.StoreSeq == 0 || req.StoreSeq > workers {
			t.Fatalf("StoreSeq = %d, want 1..%d", req.StoreSeq, workers)
		}
		if seen[req.StoreSeq] {
			t.Fatalf("duplicate StoreSeq %d", req.StoreSeq)
		}
		seen[req.StoreSeq] = true
	}
	if got := len(drainUpload(in)); got != workers {
		t.Fatalf("upload got %d, want %d", got, workers)
	}
}
