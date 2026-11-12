package manager

import (
	"context"
	"errors"
	"path/filepath"
	"testing"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func backfillScalar(num int64) *recordv1.Record {
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

func backfillHistogram(num int64) *recordv1.Record {
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

func backfillStart() *recordv1.Record {
	return &recordv1.Record{
		Num:        record.NumStart,
		RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}},
	}
}

func backfillFinish() *recordv1.Record {
	return &recordv1.Record{
		Num:        record.NumFinish,
		RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{}},
	}
}

func backfillEndSave(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{Save: &savev1.SaveRecord{
			Name:   "end.txt",
			Policy: savev1.SavePolicy_SAVE_POLICY_END,
			Type:   savev1.SaveType_SAVE_TYPE_CUSTOM,
		}},
	}
}

// newBackfillScanSession 建一个只用于扫描测试的 store。
func newBackfillScanSession(t *testing.T, records []*recordv1.Record) *Session {
	t.Helper()
	cfg := config.Defaults()
	cfg.RunID = "backfill-scan"
	cfg.RunDir = t.TempDir()
	cfg.RecordBatch = 16
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
	return &Session{cfg: cfg, store: async, runFile: path}
}

func TestScanBackfillFiltersRecordsAndRespectsFlushedBound(t *testing.T) {
	s := newBackfillScanSession(t, []*recordv1.Record{
		backfillStart(),
		backfillScalar(1),
		backfillEndSave(2),
		backfillHistogram(3),
		backfillFinish(),
	})
	_, flushedOff := s.store.Flushed()

	cursor, err := s.newBackfillCursor(2)
	if err != nil {
		t.Fatalf("newBackfillCursor: %v", err)
	}
	if cursor.seq != 1 {
		t.Fatalf("cursor seq = %d, want 1", cursor.seq)
	}
	records, _, err := s.scanBackfill(&cursor, flushedOff, 16)
	if err != nil {
		t.Fatalf("scanBackfill: %v", err)
	}
	if len(records) != 2 || records[0].GetNum() != 1 || records[1].GetNum() != 3 {
		t.Fatalf("records = %v, want scalar 1 then histogram 3", backfillNums(records))
	}
	if cursor.seq != 5 {
		t.Fatalf("cursor seq = %d, want 5 after full scan", cursor.seq)
	}

	// 上界落在最后一条记录中间：不能消费该记录，也不能越过上界。
	partial := backfillCursor{offset: store.HeaderLen}
	records, _, err = s.scanBackfill(&partial, flushedOff-1, 16)
	if err != nil {
		t.Fatalf("partial scan: %v", err)
	}
	if partial.seq != 4 {
		t.Fatalf("partial cursor seq = %d, want 4", partial.seq)
	}
	if len(records) != 2 {
		t.Fatalf("partial records = %v, want the two uploadable records", backfillNums(records))
	}
	records, _, err = s.scanBackfill(&partial, flushedOff, 16)
	if err != nil {
		t.Fatalf("resume scan: %v", err)
	}
	if len(records) != 0 || partial.seq != 5 {
		t.Fatalf("resume scan = %v / seq %d, want no uploadable records and seq 5", backfillNums(records), partial.seq)
	}
}

func commitBackfill(t *testing.T, s *Session, records ...*recordv1.Record) {
	t.Helper()
	tk, err := s.intake.Reserve(context.Background(), demandFor(records, s.cfg))
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := s.intake.Commit(tk, records); err != nil {
		t.Fatalf("Commit: %v", err)
	}
}

func TestBackfillUploadsPausedRangeInOrder(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 2, blockUpload: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	commitBackfill(t, h.session, backfillScalar(1))
	commitBackfill(t, h.session, backfillScalar(2))
	commitBackfill(t, h.session, backfillScalar(3))
	commitBackfill(t, h.session, backfillScalar(4))
	if paused, seq := h.session.intake.PauseState(); !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}

	close(h.gate)
	if err := h.session.store.Flush(ctx); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if err := h.session.runBackfill(ctx); err != nil {
		t.Fatalf("runBackfill: %v", err)
	}
	if paused, seq := h.session.intake.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState after backfill = (%v, %d), want resumed", paused, seq)
	}

	if err := h.session.intake.Close(ctx); err != nil {
		t.Fatalf("intake.Close: %v", err)
	}
	h.session.uploader.RequestFinish()
	if err := h.session.uploader.Wait(ctx); err != nil {
		t.Fatalf("uploader.Wait: %v", err)
	}
	got := h.uploads.snapshot()
	if len(got) != 4 {
		t.Fatalf("uploaded %d records, want 4", len(got))
	}
	for i, rec := range got {
		if rec.GetNum() != int64(i+1) {
			t.Fatalf("upload order = %v, want 1..4", backfillNums(got))
		}
	}
	if stored := readStore(t, h.runFile); len(stored) != 4 {
		t.Fatalf("store has %d records, want 4", len(stored))
	}
}

func TestBackfillSkipsEndSaveUntilClose(t *testing.T) {
	// 预算 2 条：END 与触发暂停的标量同批提交，延迟记录停在暂停区间。
	h := newFinishHarness(t, finishOptions{budgetRecords: 2, blockUpload: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	commitBackfill(t, h.session, backfillScalar(1))
	commitBackfill(t, h.session, backfillScalar(2))
	// 同一批里标量触发暂停，END 保存由 intake 暂存到 Close。
	commitBackfill(t, h.session, backfillScalar(3), backfillEndSave(4))
	commitBackfill(t, h.session, backfillScalar(5))
	if paused, seq := h.session.intake.PauseState(); !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}

	close(h.gate)
	if err := h.session.store.Flush(ctx); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	if err := h.session.runBackfill(ctx); err != nil {
		t.Fatalf("runBackfill: %v", err)
	}
	for _, rec := range h.uploads.snapshot() {
		if record.TypeOf(rec) == record.TypeSave {
			t.Fatal("END save must not be uploaded by backfill")
		}
	}

	if err := h.session.intake.Close(ctx); err != nil {
		t.Fatalf("intake.Close: %v", err)
	}
	h.session.uploader.RequestFinish()
	if err := h.session.uploader.Wait(ctx); err != nil {
		t.Fatalf("uploader.Wait: %v", err)
	}

	var scalars, saves int
	for _, rec := range h.uploads.snapshot() {
		switch record.TypeOf(rec) {
		case record.TypeScalar:
			scalars++
		case record.TypeSave:
			saves++
		}
	}
	if scalars != 4 || saves != 1 {
		t.Fatalf("uploaded scalars=%d saves=%d, want 4/1", scalars, saves)
	}
}

func TestFinishCompletesPausedBackfill(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 2, blockUpload: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	commitBackfill(t, h.session, backfillScalar(1))
	commitBackfill(t, h.session, backfillScalar(2))
	commitBackfill(t, h.session, backfillScalar(3))
	commitBackfill(t, h.session, backfillScalar(4))
	close(h.gate)

	if err := h.session.finish(ctx, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("finish: %v", err)
	}
	waitDrained(t, h.session)

	if paused, seq := h.session.intake.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState after finish = (%v, %d), want resumed", paused, seq)
	}
	var scalars int
	for _, rec := range h.uploads.snapshot() {
		if record.TypeOf(rec) == record.TypeScalar {
			scalars++
		}
	}
	if scalars != 4 {
		t.Fatalf("uploaded scalars = %d, want 4", scalars)
	}
	if stored := readStore(t, h.runFile); len(stored) != 5 {
		t.Fatalf("store has %d records, want 4 scalars + finish", len(stored))
	}
}

// TestFinishUploadsEndSaveAfterPausedBackfill 唯一记录预算被在途标量占满时，END 保存
// 与暂停记录共存：finish 前不上传，finish 后恰好一次。
func TestFinishUploadsEndSaveAfterPausedBackfill(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 1, blockUpload: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	commitBackfill(t, h.session, backfillScalar(1))  // 占满唯一记录预算
	commitBackfill(t, h.session, backfillEndSave(2)) // 延迟到 Close
	commitBackfill(t, h.session, backfillScalar(3))  // 预算已满 → 暂停
	if paused, seq := h.session.intake.PauseState(); !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}
	for _, rec := range h.uploads.snapshot() {
		if record.TypeOf(rec) == record.TypeSave {
			t.Fatal("END save must not be uploaded before finish")
		}
	}

	close(h.gate)
	if err := h.session.finish(ctx, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("finish: %v", err)
	}
	waitDrained(t, h.session)

	uploads := h.uploads.snapshot()
	var scalars, saves int
	for _, rec := range uploads {
		switch record.TypeOf(rec) {
		case record.TypeScalar:
			scalars++
		case record.TypeSave:
			saves++
		}
	}
	if scalars != 2 || saves != 1 {
		t.Fatalf("uploaded scalars=%d saves=%d, want 2/1", scalars, saves)
	}
	if record.TypeOf(uploads[len(uploads)-1]) != record.TypeSave {
		t.Fatalf("upload order = %v, want the END save last", backfillNums(uploads))
	}
	if paused, seq := h.session.intake.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState after finish = (%v, %d), want resumed", paused, seq)
	}
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm after finish: %v", err)
	}
}

func TestCompleteBackfillRejectsDegradedStore(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 1})
	ctx, cancel := context.WithTimeout(context.Background(), time.Second)
	defer cancel()

	// 单批两条记录超过预算，首次提交即发布暂停。
	commitBackfill(t, h.session, backfillScalar(1), backfillScalar(2))
	if paused, seq := h.session.intake.PauseState(); !paused || seq != 1 {
		t.Fatalf("PauseState = (%v, %d), want (true, 1)", paused, seq)
	}
	h.session.storeDegraded.Store(true)

	err := h.session.completeBackfill(ctx)
	if !errors.Is(err, ErrRunFailed) {
		t.Fatalf("completeBackfill err = %v, want ErrRunFailed", err)
	}
}

func backfillNums(records []*recordv1.Record) []int64 {
	nums := make([]int64, len(records))
	for i, rec := range records {
		nums[i] = rec.GetNum()
	}
	return nums
}

// 补投协程由暂停通知唤醒，恢复后自行停止。
func TestBackfillSupervisorResumesAfterPause(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 2, blockUpload: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	h.session.startBackfill()
	commitBackfill(t, h.session, backfillScalar(1))
	commitBackfill(t, h.session, backfillScalar(2))
	commitBackfill(t, h.session, backfillScalar(3))
	commitBackfill(t, h.session, backfillScalar(4))
	// 暂停期间的 END 保存不占预算，后台补投仍能推进。
	commitBackfill(t, h.session, backfillEndSave(5))
	close(h.gate)
	h.session.wakeBackfill()

	deadline := time.Now().Add(3 * time.Second)
	for {
		if paused, _ := h.session.intake.PauseState(); !paused {
			break
		}
		if time.Now().After(deadline) {
			t.Fatal("supervisor did not resume uploads")
		}
		time.Sleep(5 * time.Millisecond)
	}

	h.session.stopBackfill()
	if err := h.session.intake.Close(ctx); err != nil {
		t.Fatalf("intake.Close: %v", err)
	}
	h.session.uploader.RequestFinish()
	if err := h.session.uploader.Wait(ctx); err != nil {
		t.Fatalf("uploader.Wait: %v", err)
	}
	if got := len(h.uploads.snapshot()); got != 5 {
		t.Fatalf("uploaded %d records, want 5 (4 scalars + END save)", got)
	}
}
