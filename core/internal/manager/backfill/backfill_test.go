package backfill

import (
	"context"
	"path/filepath"
	"sync/atomic"
	"testing"
	"time"

	"go.uber.org/goleak"
	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
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
		async.In() <- store.Request{StoreSeq: uint64(i + 1), Payload: payload}
	}
	if err := async.FlushForTest(); err != nil {
		t.Fatalf("FlushForTest: %v", err)
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
	var storedSeqs []uint64
	w := New(Options{
		Store:     async,
		RunFile:   path,
		BatchSize: 16,
		OnStored: func(_ *recordv1.Record, seq uint64) {
			storedSeqs = append(storedSeqs, seq)
		},
	})
	_, flushedOff := async.Flushed()

	cursor, err := w.newCursor(2)
	if err != nil {
		t.Fatalf("newCursor: %v", err)
	}
	if cursor.Seq != 1 {
		t.Fatalf("cursor seq = %d, want 1", cursor.Seq)
	}
	records, _, err := w.scan(&cursor, flushedOff, 16)
	if err != nil {
		t.Fatalf("scan: %v", err)
	}
	if len(records) != 2 || records[0].GetNum() != 1 || records[1].GetNum() != 3 {
		t.Fatalf("records = %+v, want scalar 1 then histogram 3", records)
	}
	if cursor.Seq != 5 {
		t.Fatalf("cursor seq = %d, want 5 after full scan", cursor.Seq)
	}
	// newCursor(2) 先扫描了 seq 1，scan 再扫描 seq 2..5。
	if len(storedSeqs) != 5 {
		t.Fatalf("onStored calls = %d, want 5 (seqs 1..5)", len(storedSeqs))
	}

	// 上界落在最后一条记录中间：不能消费该记录，也不能越过上界。
	partial := Cursor{Offset: store.HeaderLen}
	records, _, err = w.scan(&partial, flushedOff-1, 16)
	if err != nil {
		t.Fatalf("partial scan: %v", err)
	}
	if partial.Seq != 4 || len(records) != 2 {
		t.Fatalf("partial = (seq=%d, len=%d), want (4, 2)", partial.Seq, len(records))
	}
	records, _, err = w.scan(&partial, flushedOff, 16)
	if err != nil {
		t.Fatalf("resume scan: %v", err)
	}
	if len(records) != 0 || partial.Seq != 5 {
		t.Fatalf("resume = (len=%d, seq=%d), want (0, 5)", len(records), partial.Seq)
	}
}

func TestCompleteRejectsDegradedOrSuspended(t *testing.T) {
	async, path := newScanStore(t, nil)
	b := budget.New(budget.Limits{Records: 1})
	intake, err := buffer.New(buffer.Options{
		Store:       async.In(),
		UploadQueue: 16,
		Budget:      b,
		Flushed: func() uint64 {
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

	degraded.Store(false)
	wSuspended := New(Options{
		Intake:      intake,
		Store:       async,
		RunFile:     path,
		IsSuspended: func() bool { return true },
	})
	if err := wSuspended.Complete(ctx); err == nil {
		t.Fatal("expected error when site is suspended")
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
