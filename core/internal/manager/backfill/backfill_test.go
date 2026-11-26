package backfill

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"go.uber.org/goleak"
	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager/progress"
	"github.com/swanhubx/swanlab/core/internal/manager/resume"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

func testScalar(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Scalar{Scalar: &datav1.ScalarRecord{
			Key:   "loss",
			Step:  num,
			Type:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
			Value: &datav1.ScalarValue{Number: float64(num)},
		}},
	}
}

func testHistogram(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Histogram{Histogram: &datav1.HistogramRecord{
			Key:   "weights",
			Step:  num,
			Type:  columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM,
			Value: &datav1.HistogramValue{Values: []uint64{1, 2}, Bins: []float64{0, 1, 2}},
		}},
	}
}

func testStart() *recordv1.Record {
	return &recordv1.Record{
		Num:        record.NumStart,
		RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}},
	}
}

func testFinish() *recordv1.Record {
	return &recordv1.Record{
		Num:        record.NumFinish,
		RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{}},
	}
}

func testEndSave(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{Save: &savev1.SaveRecord{
			Name:   "end.txt",
			Policy: savev1.SavePolicy_SAVE_POLICY_END,
			Type:   savev1.SaveType_SAVE_TYPE_CUSTOM,
		}},
	}
}

func newScanStore(t *testing.T, records []*recordv1.Record) (*store.AsyncWriter, string) {
	t.Helper()
	path := filepath.Join(t.TempDir(), "run-scan.swanlab")
	writer, err := store.Create(path)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{Queue: 64, IndexStride: 1})
	if err != nil {
		t.Fatalf("store.NewAsyncWriter: %v", err)
	}
	t.Cleanup(func() { _ = async.Close() })
	for i, rec := range records {
		payload, err := proto.Marshal(rec)
		if err != nil {
			t.Fatalf("marshal: %v", err)
		}
		async.In() <- store.Request{StoreSeq: int64(i + 1), Payload: payload}
	}
	if err := async.Flush(context.Background()); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	return async, path
}

func TestScanFiltersRecordsAndRespectsFlushedBound(t *testing.T) {
	async, path := newScanStore(t, []*recordv1.Record{
		testStart(),
		testScalar(1),
		testEndSave(2),
		testHistogram(3),
		testFinish(),
	})
	var storedSeqs []int64
	w := New(Options{
		Store:     async,
		RunFile:   path,
		BatchSize: 16,
		OnStored: func(_ *recordv1.Record, seq int64) {
			storedSeqs = append(storedSeqs, seq)
		},
	})
	_, flushedOff := async.Flushed()

	ctx := context.Background()
	cursor, err := w.newCursor(ctx, 2)
	if err != nil {
		t.Fatalf("newCursor: %v", err)
	}
	if cursor.Seq != 1 {
		t.Fatalf("cursor seq = %d, want 1", cursor.Seq)
	}
	records, _, err := w.scan(ctx, &cursor, flushedOff, 16)
	if err != nil {
		t.Fatalf("scan: %v", err)
	}
	if len(records) != 2 || records[0].GetNum() != 1 || records[1].GetNum() != 3 {
		t.Fatalf("records = %+v, want scalar 1 then histogram 3", records)
	}
	if cursor.Seq != 5 {
		t.Fatalf("cursor seq = %d, want 5 after full scan", cursor.Seq)
	}
	// newCursor(2) 仅推进帧偏移、不触发 onStored；scan 扫描 seq 2..5。
	if len(storedSeqs) != 4 || storedSeqs[0] != 2 || storedSeqs[3] != 5 {
		t.Fatalf("onStored seqs = %v, want [2 3 4 5]", storedSeqs)
	}

	// 上界落在最后一条记录中间：不能消费该记录，也不能越过上界。
	partial := Cursor{Offset: store.HeaderLen}
	records, _, err = w.scan(ctx, &partial, flushedOff-1, 16)
	if err != nil {
		t.Fatalf("partial scan: %v", err)
	}
	if partial.Seq != 4 || len(records) != 2 {
		t.Fatalf("partial = (seq=%d, len=%d), want (4, 2)", partial.Seq, len(records))
	}
	records, _, err = w.scan(ctx, &partial, flushedOff, 16)
	if err != nil {
		t.Fatalf("resume scan: %v", err)
	}
	if len(records) != 0 || partial.Seq != 5 {
		t.Fatalf("resume = (len=%d, seq=%d), want (0, 5)", len(records), partial.Seq)
	}
}

func TestCompleteRejectsDegradedStore(t *testing.T) {
	async, path := newScanStore(t, nil)
	b := budget.New(budget.Limits{Records: 1})
	intake, err := buffer.New(buffer.Options{
		Store:       async.In(),
		UploadQueue: 16,
		Budget:      b,
		Flushed: func() int64 {
			seq, _ := async.Flushed()
			return seq
		},
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}
	defer func() { _ = intake.Close(context.Background()) }()

	recs := []*recordv1.Record{testScalar(1), testScalar(2)}
	tk, err := intake.Reserve(context.Background(), buffer.Demand{Store: 2, Upload: 2, Bytes: record.Size(recs[0]) + record.Size(recs[1])})
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := intake.Commit(tk, recs); err != nil {
		t.Fatalf("Commit: %v", err)
	}
	if paused, seq := intake.PauseState(); !paused || seq != 1 {
		t.Fatalf("PauseState = (%v, %d), want (true, 1)", paused, seq)
	}

	var degraded atomic.Bool
	degraded.Store(true)
	w := New(Options{
		Intake:   intake,
		Store:    async,
		RunFile:  path,
		Degraded: &degraded,
	})
	ctx, cancel := context.WithTimeout(context.Background(), time.Second)
	defer cancel()
	if err := w.Complete(ctx); err == nil {
		t.Fatal("expected error when store is degraded")
	}
}

func TestCompleteWaitsForSiteRecovery(t *testing.T) {
	t.Run("timeout_when_suspended", func(t *testing.T) {
		async, path := newScanStore(t, []*recordv1.Record{testScalar(1), testScalar(2)})
		b := budget.New(budget.Limits{Records: 10})
		intake, err := buffer.New(buffer.Options{
			Store:       async.In(),
			UploadQueue: 16,
			Budget:      b,
			StartSeq:    2,
			Flushed: func() int64 {
				seq, _ := async.Flushed()
				return seq
			},
		})
		if err != nil {
			t.Fatalf("buffer.New: %v", err)
		}
		defer func() { _ = intake.Close(context.Background()) }()
		intake.PauseAt(1)

		w := New(Options{
			Intake:      intake,
			Store:       async,
			RunFile:     path,
			IsSuspended: func() bool { return true },
		})
		ctx, cancel := context.WithTimeout(context.Background(), 50*time.Millisecond)
		defer cancel()
		if err := w.Complete(ctx); !errors.Is(err, context.DeadlineExceeded) {
			t.Fatalf("Complete err = %v, want context.DeadlineExceeded", err)
		}
	})

	t.Run("resumes_on_health_notification", func(t *testing.T) {
		async, path := newScanStore(t, []*recordv1.Record{testScalar(1), testScalar(2), testScalar(3)})
		b := budget.New(budget.Limits{Records: 10})
		intake, err := buffer.New(buffer.Options{
			Store:       async.In(),
			UploadQueue: 16,
			Budget:      b,
			StartSeq:    3,
			Flushed: func() int64 {
				seq, _ := async.Flushed()
				return seq
			},
		})
		if err != nil {
			t.Fatalf("buffer.New: %v", err)
		}
		intake.PauseAt(1)
		if err := async.Close(); err != nil {
			t.Fatalf("async.Close: %v", err)
		}

		var suspended atomic.Bool
		suspended.Store(true)
		updates := make(chan struct{}, 1)
		var unsubscribed atomic.Bool

		w := New(Options{
			Intake:      intake,
			Store:       async,
			RunFile:     path,
			BatchSize:   16,
			IsSuspended: suspended.Load,
			SubscribeHealth: func() (<-chan struct{}, func()) {
				return updates, func() { unsubscribed.Store(true) }
			},
		})

		ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
		defer cancel()

		var uploaded []*recordv1.Record
		var mu sync.Mutex
		drainDone := make(chan struct{})
		go func() {
			defer close(drainDone)
			for rec := range intake.UploadOut() {
				mu.Lock()
				uploaded = append(uploaded, rec)
				mu.Unlock()
				b.Release(1, record.Size(rec))
			}
		}()

		done := make(chan error, 1)
		go func() { done <- w.Complete(ctx) }()

		select {
		case err := <-done:
			t.Fatalf("Complete returned early while site suspended: %v", err)
		case <-time.After(30 * time.Millisecond):
		}

		suspended.Store(false)
		updates <- struct{}{}

		if err := <-done; err != nil {
			t.Fatalf("Complete after site recovery: %v", err)
		}
		if err := intake.Close(ctx); err != nil {
			t.Fatalf("intake.Close: %v", err)
		}
		<-drainDone

		if !unsubscribed.Load() {
			t.Fatal("SubscribeHealth unsubscribe was not called")
		}
		mu.Lock()
		defer mu.Unlock()
		if len(uploaded) != 3 {
			t.Fatalf("uploaded %d records, want 3", len(uploaded))
		}
		for i, rec := range uploaded {
			if rec.GetNum() != int64(i+1) {
				t.Fatalf("uploaded[%d].Num = %d, want %d", i, rec.GetNum(), i+1)
			}
		}
	})
}

func TestCompleteRetriesOnOverloadedMidBackfill(t *testing.T) {
	const total = 8
	recs := make([]*recordv1.Record, 0, total)
	for i := 1; i <= total; i++ {
		recs = append(recs, testScalar(int64(i)))
	}
	async, path := newScanStore(t, recs)
	b := budget.New(budget.Limits{})
	intake, err := buffer.New(buffer.Options{
		Store:       async.In(),
		UploadQueue: 2,
		QueueWait:   10 * time.Millisecond,
		Budget:      b,
		StartSeq:    total,
		Flushed: func() int64 {
			seq, _ := async.Flushed()
			return seq
		},
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}
	intake.PauseAt(1)
	if err := async.Close(); err != nil {
		t.Fatalf("async.Close: %v", err)
	}

	w := New(Options{
		Intake:    intake,
		Store:     async,
		RunFile:   path,
		BatchSize: 2,
	})

	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	var uploaded []*recordv1.Record
	drainDone := make(chan struct{})
	go func() {
		defer close(drainDone)
		first := true
		for rec := range intake.UploadOut() {
			uploaded = append(uploaded, rec)
			if first {
				first = false
				// 暂停消费超过 QueueWait，使 Replay 触发 ErrOverloaded 并从当前游标重试。
				time.Sleep(50 * time.Millisecond)
			}
		}
	}()

	if err := w.Complete(ctx); err != nil {
		t.Fatalf("Complete with mid-backfill overload: %v", err)
	}
	if err := intake.Close(ctx); err != nil {
		t.Fatalf("intake.Close: %v", err)
	}
	<-drainDone

	if len(uploaded) != total {
		t.Fatalf("uploaded %d records, want %d", len(uploaded), total)
	}
	for i, rec := range uploaded {
		if rec.GetNum() != int64(i+1) {
			t.Fatalf("uploaded[%d].Num = %d, want %d", i, rec.GetNum(), i+1)
		}
	}
}

func TestNilWorkerIsSafe(t *testing.T) {
	var w *Worker
	w.Start()
	w.Wake()
	w.CancelAttempt()
	w.Stop()
	if err := w.RunOnce(context.Background()); err != nil {
		t.Fatalf("RunOnce on nil worker = %v", err)
	}
	if err := w.Complete(context.Background()); err != nil {
		t.Fatalf("Complete on nil worker = %v", err)
	}
}

func TestPausedWritesKeepProgressTrackerBoundedAndResumeWatermark(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run-paused-progress.swanlab")
	writer, err := store.Create(path)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{Queue: 512, IndexStride: 256})
	if err != nil {
		t.Fatalf("store.NewAsyncWriter: %v", err)
	}
	defer func() { _ = async.Close() }()

	tracker := progress.NewTracker(0)
	b := budget.New(budget.Limits{Records: 2})
	intake, err := buffer.New(buffer.Options{
		Store:       async.In(),
		UploadQueue: 64,
		Budget:      b,
		Flushed: func() int64 {
			seq, _ := async.Flushed()
			return seq
		},
		OnStored: tracker.OnStored,
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()

	// 先写入并确认 2 条实时记录，保持其占用 2 条预算。
	live := []*recordv1.Record{testScalar(1), testScalar(2)}
	tk, err := intake.Reserve(ctx, buffer.Demand{
		Store:  2,
		Upload: 2,
		Bytes:  record.Size(live[0]) + record.Size(live[1]),
	})
	if err != nil {
		t.Fatalf("Reserve live: %v", err)
	}
	if err = intake.Commit(tk, live); err != nil {
		t.Fatalf("Commit live: %v", err)
	}

	// 确认前 2 条记录上传完成（暂不释放预算，让后续记录触发暂停）。
	for i := 0; i < 2; i++ {
		rec := <-intake.UploadOut()
		tracker.Resolve([]*recordv1.Record{rec})
	}
	if got := tracker.Watermark(); got != 2 {
		t.Fatalf("live watermark = %d, want 2", got)
	}

	// 预算仍占满，写入 10_000 条暂停记录与 FinishRecord。
	const pausedScalars = 10_000
	for i := 0; i < pausedScalars; i += 200 {
		batch := make([]*recordv1.Record, 0, 200)
		var bytes int64
		for j := 0; j < 200; j++ {
			rec := testScalar(int64(3 + i + j))
			batch = append(batch, rec)
			bytes += record.Size(rec)
		}
		var ptk buffer.Ticket
		ptk, err = intake.Reserve(ctx, buffer.Demand{Store: len(batch), Upload: len(batch), Bytes: bytes})
		if err != nil {
			t.Fatalf("Reserve paused batch: %v", err)
		}
		if err = intake.Commit(ptk, batch); err != nil {
			t.Fatalf("Commit paused batch: %v", err)
		}
	}
	ftk, err := intake.Reserve(ctx, buffer.Demand{Store: 1})
	if err != nil {
		t.Fatalf("Reserve finish: %v", err)
	}
	if err = intake.Commit(ftk, []*recordv1.Record{testFinish()}); err != nil {
		t.Fatalf("Commit finish: %v", err)
	}

	if paused, seq := intake.PauseState(); !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}
	// FinishRecord 已写入，但在补投确认前水位必须停在 2。
	if got := tracker.Watermark(); got != 2 {
		t.Fatalf("watermark before backfill = %d, want 2", got)
	}

	// 释放前面 2 条记录占用的预算，模拟消费上传通道并执行补投。
	b.Release(2, record.Size(live[0])+record.Size(live[1]))
	if err := async.Close(); err != nil {
		t.Fatalf("async.Close: %v", err)
	}

	drainDone := make(chan struct{})
	go func() {
		defer close(drainDone)
		for rec := range intake.UploadOut() {
			tracker.Resolve([]*recordv1.Record{rec})
			b.Release(1, record.Size(rec))
		}
	}()

	w := New(Options{
		Intake:    intake,
		Store:     async,
		RunFile:   path,
		BatchSize: 64,
		OnStored:  tracker.OnStored,
	})
	if err := w.Complete(ctx); err != nil {
		t.Fatalf("Complete: %v", err)
	}
	if err := intake.Close(ctx); err != nil {
		t.Fatalf("intake.Close: %v", err)
	}
	<-drainDone

	wantWatermark := int64(2 + pausedScalars + 1)
	if got := tracker.Watermark(); got != wantWatermark {
		t.Fatalf("final watermark = %d, want %d", got, wantWatermark)
	}
}

func TestNewCursorPreservesPrePauseProgressAndAdvancesCheckpoint(t *testing.T) {
	pre1, pre2 := testScalar(1), testScalar(2)
	post3, post4 := testScalar(3), testScalar(4)
	async, path := newScanStore(t, []*recordv1.Record{pre1, pre2, post3, post4})
	_, flushedOff := async.Flushed()

	tracker := progress.NewTracker(0)
	// 暂停前记录仍停留在上传管线中，尚未确认。
	tracker.OnStored(pre1, 1)
	tracker.OnStored(pre2, 2)

	w := New(Options{
		Store:     async,
		RunFile:   path,
		BatchSize: 16,
		OnStored:  tracker.OnStored,
	})
	ctx := context.Background()
	cursor, err := w.newCursor(ctx, 3)
	if err != nil {
		t.Fatalf("newCursor: %v", err)
	}
	replayed, _, err := w.scan(ctx, &cursor, flushedOff, 16)
	if err != nil {
		t.Fatalf("scan: %v", err)
	}
	if len(replayed) != 2 {
		t.Fatalf("replayed len = %d, want 2", len(replayed))
	}

	// 先确认暂停后的补投记录，再确认暂停前的原始指针：水位必须追平到 4。
	tracker.Resolve(replayed)
	if got := tracker.Watermark(); got != 0 {
		t.Fatalf("watermark before pre-pause confirm = %d, want 0", got)
	}
	tracker.Resolve([]*recordv1.Record{pre1, pre2})
	if got := tracker.Watermark(); got != 4 {
		t.Fatalf("final watermark = %d, want 4", got)
	}

	cpOff, err := store.EndOffset(path, async.Index(), tracker.Watermark(), flushedOff)
	if err != nil {
		t.Fatalf("EndOffset: %v", err)
	}
	cpFile := filepath.Join(filepath.Dir(path), "checkpoint.json")
	if err = store.WriteCheckpointFast(cpFile, store.Checkpoint{
		RunID:    "run-cp",
		StoreSeq: tracker.Watermark(),
		Offset:   cpOff,
	}); err != nil {
		t.Fatalf("WriteCheckpointFast: %v", err)
	}
	cp, ok, err := store.LoadCheckpoint(cpFile, "run-cp")
	if err != nil || !ok || cp.StoreSeq != 4 || cp.Offset != flushedOff {
		t.Fatalf("checkpoint = (%+v, %v, %v), want StoreSeq=4 Offset=%d", cp, ok, err, flushedOff)
	}
}

func TestScanAndSkipRespectContextCancellation(t *testing.T) {
	async, path := newScanStore(t, []*recordv1.Record{testStart(), testScalar(1), testScalar(2), testScalar(3)})
	_, flushedOff := async.Flushed()
	w := New(Options{Store: async, RunFile: path, BatchSize: 16})

	canceled, cancel := context.WithCancel(context.Background())
	cancel()

	if _, err := w.newCursor(canceled, 4); !errors.Is(err, context.Canceled) {
		t.Fatalf("newCursor canceled err = %v, want context.Canceled", err)
	}
	c := Cursor{Offset: store.HeaderLen}
	if _, _, err := w.scan(canceled, &c, flushedOff, 16); !errors.Is(err, context.Canceled) {
		t.Fatalf("scan canceled err = %v, want context.Canceled", err)
	}
}

func TestResume100kRecordsPauseNewCursorDecodesZero(t *testing.T) {
	const total = 100_000
	cfg := config.Defaults()
	cfg.RunID = "run-100k"
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	runFile, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}
	if err = os.MkdirAll(filepath.Dir(runFile), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	writer, err := store.Create(runFile)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	startPayload, err := proto.Marshal(testStart())
	if err != nil {
		_ = writer.Close()
		t.Fatalf("marshal start: %v", err)
	}
	if err = writer.Write(startPayload); err != nil {
		_ = writer.Close()
		t.Fatalf("write start: %v", err)
	}
	scalarPayload, err := proto.Marshal(testScalar(1))
	if err != nil {
		_ = writer.Close()
		t.Fatalf("marshal scalar: %v", err)
	}
	for i := 2; i <= total; i++ {
		if err = writer.Write(scalarPayload); err != nil {
			_ = writer.Close()
			t.Fatalf("write[%d]: %v", i, err)
		}
	}
	if err = writer.Close(); err != nil {
		t.Fatalf("writer.Close: %v", err)
	}

	local, err := metrics.Rebuild(runFile)
	if err != nil {
		t.Fatalf("Rebuild: %v", err)
	}
	opened, err := resume.OpenStore(resume.StoreOptions{Run: run, Resume: local, Queue: 16})
	if err != nil {
		t.Fatalf("OpenStore: %v", err)
	}
	defer func() { _ = opened.Writer.Close() }()

	// 续传后立即写入第 100,001 条记录并触发暂停，此时 flushedSeq 推进至 100,001。
	opened.Writer.In() <- store.Request{StoreSeq: total + 1, Payload: scalarPayload}
	if err = opened.Writer.Flush(context.Background()); err != nil {
		t.Fatalf("Flush: %v", err)
	}

	var decodedCalls int
	w := New(Options{
		Store:   opened.Writer,
		RunFile: runFile,
		OnStored: func(*recordv1.Record, int64) {
			decodedCalls++
		},
	})
	cursor, err := w.newCursor(context.Background(), total+1)
	if err != nil {
		t.Fatalf("newCursor: %v", err)
	}
	if cursor.Seq != total || cursor.Offset != local.EndOffset {
		t.Fatalf("cursor = %+v, want Seq=%d Offset=%d", cursor, total, local.EndOffset)
	}
	if decodedCalls != 0 {
		t.Fatalf("newCursor decoded/onStored calls = %d, want 0", decodedCalls)
	}
}

func TestNewCursorDetectsCorruptPayload(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run-corrupt.swanlab")
	writer, err := store.Create(path)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	// IndexStride 设为 3，保证定位第 5 条时命中第 3 条锚点并经 skipRecords 校验损坏的第 3 条。
	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{Queue: 32, IndexStride: 3})
	if err != nil {
		_ = writer.Close()
		t.Fatalf("store.NewAsyncWriter: %v", err)
	}
	defer func() { _ = async.Close() }()

	recs := []*recordv1.Record{testStart(), testScalar(1), testScalar(2), testScalar(3), testScalar(4)}
	for i, rec := range recs {
		payload, merr := proto.Marshal(rec)
		if merr != nil {
			t.Fatalf("marshal: %v", merr)
		}
		async.In() <- store.Request{StoreSeq: int64(i + 1), Payload: payload}
	}
	if ferr := async.Flush(context.Background()); ferr != nil {
		t.Fatalf("Flush: %v", ferr)
	}

	// 找到第 3 条记录起点并原地翻转其 payload 首字节。
	reader, err := store.Open(path)
	if err != nil {
		t.Fatalf("store.Open: %v", err)
	}
	for range 2 {
		if serr := reader.Skip(); serr != nil {
			_ = reader.Close()
			t.Fatalf("Skip prefix: %v", serr)
		}
	}
	rec3PayloadOff := reader.Position() + store.HeaderLen
	_ = reader.Close()

	fp, err := os.OpenFile(path, os.O_RDWR, 0)
	if err != nil {
		t.Fatalf("OpenFile: %v", err)
	}
	var b [1]byte
	if _, rerr := fp.ReadAt(b[:], rec3PayloadOff); rerr != nil {
		_ = fp.Close()
		t.Fatalf("ReadAt: %v", rerr)
	}
	b[0] ^= 0xFF
	if _, werr := fp.WriteAt(b[:], rec3PayloadOff); werr != nil {
		_ = fp.Close()
		t.Fatalf("WriteAt: %v", werr)
	}
	_ = fp.Close()

	w := New(Options{Store: async, RunFile: path, BatchSize: 16})
	_, cursorErr := w.newCursor(context.Background(), 5)
	if !errors.Is(cursorErr, errCorrupt) || !errors.Is(cursorErr, store.ErrCorrupt) {
		t.Fatalf("newCursor err = %v, want errCorrupt wrapping store.ErrCorrupt", cursorErr)
	}
}

func TestWorkerStopBeforeStartAndConcurrentLifecycle(t *testing.T) {
	// Stop -> Start -> Stop：先 Stop 后，再 Start 不得启动后台协程或触发 onStored。
	for range 16 {
		async, path := newScanStore(t, []*recordv1.Record{testStart(), testScalar(1)})
		b := budget.New(budget.Limits{Records: 16})
		in, err := buffer.New(buffer.Options{
			Store:       async.In(),
			UploadQueue: 8,
			Budget:      b,
			Flushed: func() int64 {
				seq, _ := async.Flushed()
				return seq
			},
		})
		if err != nil {
			t.Fatalf("buffer.New: %v", err)
		}
		var storedCalls atomic.Int32
		w := New(Options{
			Intake:   in,
			Store:    async,
			RunFile:  path,
			OnStored: func(*recordv1.Record, int64) { storedCalls.Add(1) },
		})

		w.Stop()
		w.Start()
		w.Stop()
		w.Stop()
		if w.Started() {
			t.Fatal("worker must not report Started after Stop -> Start")
		}
		if storedCalls.Load() != 0 {
			t.Fatalf("storedCalls = %d, want 0 after Stop -> Start", storedCalls.Load())
		}
		_ = in.Close(context.Background())
		_ = async.Close()
	}

	// 并发 Start / Stop 在超时预算内稳定退出，且可重复 Stop、nil/未启用安全调用。
	for range 24 {
		async, path := newScanStore(t, []*recordv1.Record{testStart(), testScalar(1)})
		in, err := buffer.New(buffer.Options{
			Store:       async.In(),
			UploadQueue: 8,
			Budget:      budget.New(budget.Limits{Records: 16}),
			Flushed: func() int64 {
				seq, _ := async.Flushed()
				return seq
			},
		})
		if err != nil {
			t.Fatalf("buffer.New: %v", err)
		}
		w := New(Options{Intake: in, Store: async, RunFile: path})

		var wg sync.WaitGroup
		wg.Add(4)
		go func() { defer wg.Done(); w.Start() }()
		go func() { defer wg.Done(); w.Start() }()
		go func() { defer wg.Done(); w.Stop() }()
		go func() { defer wg.Done(); w.Stop() }()

		done := make(chan struct{})
		go func() {
			wg.Wait()
			close(done)
		}()
		select {
		case <-done:
		case <-time.After(2 * time.Second):
			t.Fatal("concurrent Worker Start/Stop timed out")
		}
		w.Stop()
		_ = in.Close(context.Background())
		_ = async.Close()
	}

	// nil 与未启用 Worker 安全空调用。
	var nilWorker *Worker
	nilWorker.Start()
	nilWorker.Stop()
	disabled := New(Options{})
	disabled.Start()
	disabled.Stop()
}
