package manager

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/cron/health"
	"github.com/swanhubx/swanlab/core/internal/manager/checkpoint"
	"github.com/swanhubx/swanlab/core/internal/manager/planner"
	"github.com/swanhubx/swanlab/core/internal/manager/progress"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func waitDrained(t *testing.T, s *Session) {
	t.Helper()
	select {
	case <-s.drained:
	case <-time.After(3 * time.Second):
		t.Fatal("session did not drain")
	}
}

func TestConfirmSucceedsAfterPausedBackfill(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 1})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	// 单批两条记录超过预算，Intake 发布暂停且两条记录只落本地。
	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("finish/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("finish/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	if paused, seq := h.session.intake.PauseState(); !paused || seq == 0 {
		t.Fatalf("pause state = %v/%d, want paused with seq", paused, seq)
	}

	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm after backfill: %v", err)
	}
	if paused, seq := h.session.intake.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState = (%v, %d), want resumed", paused, seq)
	}
}

func TestDeliverRunFinishFlushesAndIsIdempotent(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("finish/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}

	rec := finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")
	if err := h.mgr.DeliverRunFinish(ctx, h.id, rec); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	// 重复投递必须幂等，且不能产生第二条 FinishRecord。
	if err := h.mgr.DeliverRunFinish(ctx, h.id, rec); err != nil {
		t.Fatalf("repeated DeliverRunFinish: %v", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining && state != RunStateDrained {
		t.Fatalf("state after DeliverRunFinish = %s, want draining or drained", state)
	}
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	waitDrained(t, h.session)

	records := readStore(t, h.runFile)
	if len(records) == 0 {
		t.Fatal("run file is empty")
	}
	last := records[len(records)-1]
	if last.GetNum() != record.NumFinish {
		t.Fatalf("last record num = %d, want %d", last.GetNum(), record.NumFinish)
	}
	if last.GetFinish().GetState() != runv1.RunState_RUN_STATE_FINISHED {
		t.Fatalf("finish state = %v, want FINISHED", last.GetFinish().GetState())
	}
	sawFinish := 0
	for _, item := range records {
		if item.GetNum() == record.NumFinish {
			sawFinish++
		}
	}
	if sawFinish != 1 {
		t.Fatalf("stored %d finish records, want 1", sawFinish)
	}

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("finish/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("upsert after confirm err=%v, want ErrSessionNotFound", err)
	}
}

func TestDeliverRunFinishUploadsErrorLogBeforeDrain(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	rec := finishRecord(runv1.RunState_RUN_STATE_CRASHED, "boom")
	if err := h.mgr.DeliverRunFinish(ctx, h.id, rec); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	waitDrained(t, h.session)

	uploads := h.uploads.snapshot()
	if len(uploads) != 1 || uploads[0].GetLog() == nil {
		t.Fatalf("uploads = %v, want one error log", uploads)
	}
	log := uploads[0].GetLog()
	if log.GetLine() != "boom" {
		t.Fatalf("error log line = %q, want boom", log.GetLine())
	}
	if log.GetLevel() != terminalv1.LogLevel_LOG_LEVEL_ERROR {
		t.Fatalf("error log level = %v, want ERROR", log.GetLevel())
	}
	if got := log.GetTimestamp().AsTime(); !got.Equal(rec.GetFinishedAt().AsTime()) {
		t.Fatalf("error log timestamp = %v, want %v", got, rec.GetFinishedAt().AsTime())
	}

	records := readStore(t, h.runFile)
	if len(records) != 1 || records[0].GetFinish().GetState() != runv1.RunState_RUN_STATE_CRASHED {
		t.Fatalf("stored records = %v, want one CRASHED finish record", records)
	}
}

func TestDeliverRunFinishUsesDefaultErrorText(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_ABORTED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	waitDrained(t, h.session)

	uploads := h.uploads.snapshot()
	if len(uploads) != 1 || uploads[0].GetLog().GetLine() != finishUnknownError {
		t.Fatalf("error log = %v, want default message %q", uploads, finishUnknownError)
	}
}

func TestDeliverRunFinishStoreFailureDegradesOnly(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	// 模拟磁盘持续写入失败：AsyncWriter 只能降级，收尾仍须在 Confirm 完成上传侧排空。
	if err := h.writer.Close(); err != nil {
		t.Fatalf("close writer: %v", err)
	}
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish with broken store: %v", err)
	}
	if !h.session.storeDegraded.Load() {
		t.Fatal("store failure must set storeDegraded")
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining && state != RunStateDrained {
		t.Fatalf("state = %s, want draining or drained", state)
	}
	// store 降级时 finishCommitted 同样置位，Confirm 仍可完成排空与释放。
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateDrained {
		t.Fatalf("state = %s, want drained", state)
	}
	if _, err := h.mgr.State(h.id); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

func TestDeliverRunFinishLocalStoreDegradedFails(t *testing.T) {
	for _, mode := range []config.Mode{config.ModeLocal, config.ModeOffline} {
		t.Run(string(mode), func(t *testing.T) {
			h := newFinishHarness(t, finishOptions{mode: mode})
			ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
			defer cancel()

			h.session.storeDegraded.Store(true)
			err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
			if !errors.Is(err, ErrRunFailed) {
				t.Fatalf("DeliverRunFinish err=%v, want ErrRunFailed", err)
			}
			if state := h.session.stateSnapshot(); state != RunStateFailed {
				t.Fatalf("state = %s, want failed", state)
			}
			if h.session.finishCommitted {
				t.Fatal("degraded store must not mark the finish record committed")
			}
			if _, err = h.mgr.Confirm(ctx, h.id); !errors.Is(err, ErrRunFailed) {
				t.Fatalf("Confirm err=%v, want ErrRunFailed", err)
			}
			if _, err = h.mgr.State(h.id); err != nil {
				t.Fatalf("session must stay registered for cleanup/retry: %v", err)
			}
		})
	}
}

func TestDeliverRunFinishLocalStoreCloseFailureFails(t *testing.T) {
	for _, mode := range []config.Mode{config.ModeLocal, config.ModeOffline} {
		t.Run(string(mode), func(t *testing.T) {
			h := newFinishHarness(t, finishOptions{mode: mode})
			ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
			defer cancel()

			// 关闭底层文件，让 AsyncWriter 在 Close 刷盘时报错并置 storeDegraded。
			if err := h.writer.Close(); err != nil {
				t.Fatalf("close writer: %v", err)
			}
			err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
			if !errors.Is(err, ErrRunFailed) {
				t.Fatalf("DeliverRunFinish err=%v, want ErrRunFailed", err)
			}
			if !h.session.storeDegraded.Load() {
				t.Fatal("store close failure must set storeDegraded")
			}
			if state := h.session.stateSnapshot(); state != RunStateFailed {
				t.Fatalf("state = %s, want failed", state)
			}
			if _, err = h.mgr.Confirm(ctx, h.id); !errors.Is(err, ErrRunFailed) {
				t.Fatalf("Confirm err=%v, want ErrRunFailed", err)
			}
		})
	}
}

func TestDeliverRunFinishLocalModesPersistAndConfirm(t *testing.T) {
	for _, mode := range []config.Mode{config.ModeLocal, config.ModeOffline} {
		t.Run(string(mode), func(t *testing.T) {
			h := newFinishHarness(t, finishOptions{mode: mode})
			ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
			defer cancel()

			if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
				scalarRecord("finish/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
			}); err != nil {
				t.Fatalf("UpsertScalars: %v", err)
			}
			if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
				t.Fatalf("DeliverRunFinish: %v", err)
			}
			waitDrained(t, h.session)

			records := readStore(t, h.runFile)
			if len(records) == 0 {
				t.Fatal("run file is empty")
			}
			if last := records[len(records)-1]; last.GetNum() != record.NumFinish {
				t.Fatalf("last record num = %d, want %d", last.GetNum(), record.NumFinish)
			}

			result, err := h.mgr.Confirm(ctx, h.id)
			if err != nil || result.ReleaseFailed {
				t.Fatalf("Confirm result=%+v err=%v, want success", result, err)
			}
			if _, err = h.mgr.State(h.id); !errors.Is(err, ErrSessionNotFound) {
				t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
			}
		})
	}
}

func TestDeliverRunFinishIntakeCloseFailureMarksFailed(t *testing.T) {
	h := newFinishHarness(t, finishOptions{
		skipStore:   true,
		uploadQueue: 2,
		queueWait:   20 * time.Millisecond,
		blockUpload: true,
	})
	setupCtx := context.Background()
	if err := h.mgr.UpsertScalars(setupCtx, h.id, []*datav1.ScalarRecord{
		scalarRecord("finish/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	select {
	case <-h.uploadEntered:
	case <-time.After(3 * time.Second):
		t.Fatal("uploader did not enter the sender")
	}
	for step := 1; step <= 16; step++ {
		err := h.mgr.UpsertScalars(setupCtx, h.id, []*datav1.ScalarRecord{
			scalarRecord("finish/loss", int64(step), columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		})
		if errors.Is(err, ErrOverloaded) {
			break
		}
		if err != nil {
			t.Fatalf("fill upload queue: %v", err)
		}
	}
	endPath := filepath.Join(t.TempDir(), "end.txt")
	if err := os.WriteFile(endPath, []byte("end"), 0o644); err != nil {
		t.Fatalf("write end file: %v", err)
	}
	if err := h.mgr.UpsertSaves(setupCtx, h.id, []*savev1.SaveRecord{{
		Name:       "end.txt",
		SourcePath: endPath,
		Policy:     savev1.SavePolicy_SAVE_POLICY_END,
		Type:       savev1.SaveType_SAVE_TYPE_CUSTOM,
	}}); err != nil {
		t.Fatalf("UpsertSaves: %v", err)
	}

	if err := h.mgr.DeliverRunFinish(setupCtx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	// 上传通道持续满时 Confirm 短超时返回 DeadlineExceeded，且保留 DRAINING 状态供重试。
	ctx, cancel := context.WithTimeout(context.Background(), 50*time.Millisecond)
	defer cancel()
	_, err := h.mgr.Confirm(ctx, h.id)
	if !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("Confirm err=%v, want context.DeadlineExceeded", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining {
		t.Fatalf("state = %s, want draining", state)
	}
}

// newBlockedStoreHarness 用容量为 1 的受控本地通道装配会话：预填一条占位后
// FinishRecord 的 Reserve 必然失败，释放占位后重试才能成功。
func newBlockedStoreHarness(t *testing.T) (*Manager, *Session, string, chan store.Request) {
	t.Helper()
	cfg := config.Defaults()
	cfg.RunID = "blocked-store"
	cfg.RecordBatch = 8
	cfg.RecordInterval = 0.01

	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.counter = counter.New(0)
	s.epoch = counter.New(0)
	s.watcher = watcher.NullWatcher{}
	s.tracker = uploader.NewTracker()

	storeCh := make(chan store.Request, 1)
	// 占满唯一槽位，模拟本地写入队列压力。
	storeCh <- store.Request{StoreSeq: 0}

	b := budget.New(budget.Limits{})
	intake, err := buffer.New(buffer.Options{
		Store:       storeCh,
		UploadQueue: 8,
		Budget:      b,
		Flushed:     func() uint64 { return 0 },
		OnAdmit:     s.tracker.AddTotalRecords,
		QueueWait:   20 * time.Millisecond,
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}
	up, err := uploader.New(uploader.Options{
		In:       intake.UploadOut(),
		Budget:   b,
		Senders:  map[record.Type]uploader.Sender{record.TypeLog: &collectingSender{}},
		Batch:    4,
		Interval: 10 * time.Millisecond,
		Tracker:  s.tracker,
		Context:  s.ctx,
		OnDrained: func() {
			s.markDrained()
		},
	})
	if err != nil {
		t.Fatalf("uploader.New: %v", err)
	}
	s.intake = intake
	s.uploader = up
	s.drainFn = up.Wait

	mgr := newTestManager()
	id, err := mgr.install(s)
	if err != nil {
		t.Fatalf("install: %v", err)
	}
	t.Cleanup(func() {
		s.cancel()
		_ = up.Close()
	})
	return mgr, s, id, storeCh
}

// FinishRecord 落盘失败后 Confirm 必须重试，不能直接上报终态并释放会话。
func TestConfirmRetriesFinishRecordCommit(t *testing.T) {
	mgr, s, id, storeCh := newBlockedStoreHarness(t)
	var releases int
	s.releaseFn = func(context.Context) error {
		releases++
		return nil
	}

	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	err := mgr.DeliverRunFinish(ctx, id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
	if !errors.Is(err, ErrRunFailed) {
		t.Fatalf("DeliverRunFinish err=%v, want ErrRunFailed", err)
	}
	if !errors.Is(err, errFinishRecord) {
		t.Fatalf("DeliverRunFinish err=%v, want finish-record failure", err)
	}
	mustState(t, mgr, id, RunStateFailed)
	if s.finishCommitted {
		t.Fatal("finish record must stay uncommitted after a failed commit")
	}
	if releases != 0 {
		t.Fatal("release hook must not run before the finish record is persisted")
	}
	if got := len(storeCh); got != 1 {
		t.Fatalf("store queue length = %d, want only the placeholder", got)
	}

	<-storeCh // 腾出唯一槽位，模拟本地队列压力消失

	result, err := mgr.Confirm(ctx, id)
	if err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if result.ReleaseFailed {
		t.Fatalf("Confirm = %+v, want success", result)
	}
	if releases != 1 {
		t.Fatalf("release hook called %d times, want 1", releases)
	}
	if got := len(storeCh); got != 1 {
		t.Fatalf("store queue length after confirm = %d, want 1 finish record", got)
	}
	req := <-storeCh
	stored := &recordv1.Record{}
	if err = proto.Unmarshal(req.Payload, stored); err != nil {
		t.Fatalf("unmarshal finish record: %v", err)
	}
	if stored.GetNum() != record.NumFinish ||
		stored.GetFinish().GetState() != runv1.RunState_RUN_STATE_FINISHED {
		t.Fatalf("stored record = %v, want FINISHED finish record", stored)
	}
	if _, err = mgr.State(id); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

// 重复投递必须从失败点续跑，且不能写入第二条 FinishRecord。
func TestDeliverRunFinishRetriesWithoutDuplicateRecord(t *testing.T) {
	mgr, s, id, storeCh := newBlockedStoreHarness(t)

	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	rec := finishRecord(runv1.RunState_RUN_STATE_CRASHED, "boom")
	if err := mgr.DeliverRunFinish(ctx, id, rec); !errors.Is(err, ErrRunFailed) {
		t.Fatalf("first DeliverRunFinish err=%v, want ErrRunFailed", err)
	}
	<-storeCh // 腾出唯一槽位

	if err := mgr.DeliverRunFinish(ctx, id, rec); err != nil {
		t.Fatalf("retry DeliverRunFinish: %v", err)
	}
	if !s.finishCommitted {
		t.Fatal("retry must commit the finish record")
	}
	if !s.hasPendingLog() {
		t.Fatal("crash error log must stay pending until Confirm")
	}
	if got := len(storeCh); got != 1 {
		t.Fatalf("store queue length = %d, want 1", got)
	}
	req := <-storeCh
	stored := &recordv1.Record{}
	if err := proto.Unmarshal(req.Payload, stored); err != nil {
		t.Fatalf("unmarshal finish record: %v", err)
	}
	if stored.GetNum() != record.NumFinish {
		t.Fatalf("stored record num = %d, want %d", stored.GetNum(), record.NumFinish)
	}

	// 收尾已完成后再投递必须幂等，不能再写一条 -2。
	if err := mgr.DeliverRunFinish(ctx, id, rec); err != nil {
		t.Fatalf("repeated DeliverRunFinish: %v", err)
	}
	select {
	case <-storeCh:
		t.Fatal("repeated DeliverRunFinish wrote a second finish record")
	default:
	}

	result, err := mgr.Confirm(ctx, id)
	if err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if result.ReleaseFailed {
		t.Fatalf("Confirm = %+v, want success", result)
	}
	if s.hasPendingLog() {
		t.Fatal("Confirm must queue and clear the crash error log")
	}
	if _, err = mgr.State(id); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

// 上传通道暂时满时崩溃 error_log 不能被丢弃：保留到 Confirm 恢复容量后补传。
func TestDeliverRunFinishRetriesCrashLogAfterCapacityFrees(t *testing.T) {
	h := newFinishHarness(t, finishOptions{
		// 有本地副本时上传通道满会直接转暂停落盘，不会返回 ErrOverloaded；
		// 这里用 skipStore 验证无本地副本时的容量阻塞与崩溃 error_log 重试。
		skipStore:   true,
		uploadQueue: 2,
		queueWait:   200 * time.Millisecond,
		blockUpload: true,
	})
	ctx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()

	// 第一条指标（含自动列）把 uploader 扣在 sender 上，后续记录只能排在通道里。
	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("crash/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	select {
	case <-h.uploadEntered:
	case <-time.After(3 * time.Second):
		t.Fatal("uploader did not enter the sender")
	}

	// uploader 阻塞期间不会再消费，持续写入直到上传通道占满。
	saturated := false
	for step := 1; step <= 16; step++ {
		err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
			scalarRecord("crash/loss", int64(step), columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		})
		if err == nil {
			continue
		}
		if !errors.Is(err, ErrOverloaded) {
			t.Fatalf("fill upload queue: %v", err)
		}
		saturated = true
		break
	}
	if !saturated {
		t.Fatal("upload queue was never saturated")
	}

	// DeliverRunFinish 只做本地收尾，不受上传通道满影响，应立即成功。
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_CRASHED, "boom")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining {
		t.Fatalf("state = %s, want draining", state)
	}

	// 上传通道仍满时 Confirm 短超时返回 DeadlineExceeded，但保留 pendingLog 供下次重试。
	shortCtx, shortCancel := context.WithTimeout(ctx, 100*time.Millisecond)
	_, err := h.mgr.Confirm(shortCtx, h.id)
	shortCancel()
	if !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("first Confirm err=%v, want context.DeadlineExceeded", err)
	}
	if !h.session.hasPendingLog() {
		t.Fatal("crash error log must be kept for retry")
	}
	for _, rec := range h.uploads.snapshot() {
		if rec.GetLog() != nil {
			t.Fatal("crash error log must not be uploaded before retry")
		}
	}

	close(h.gate) // 放行上传，让通道重新腾出槽位

	result, err := h.mgr.Confirm(ctx, h.id)
	if err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if result.ReleaseFailed {
		t.Fatalf("Confirm = %+v, want success", result)
	}
	if h.session.hasPendingLog() {
		t.Fatal("crash error log must be cleared after a successful upload")
	}
	var logs []*recordv1.Record
	for _, rec := range h.uploads.snapshot() {
		if rec.GetLog() != nil {
			logs = append(logs, rec)
		}
	}
	if len(logs) != 1 {
		t.Fatalf("uploaded %d error log(s), want 1", len(logs))
	}
	if got := logs[0].GetLog().GetLine(); got != "boom" {
		t.Fatalf("error log line = %q, want boom", got)
	}
	if _, err = h.mgr.State(h.id); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

// 上传通道持续满超过 QueueWait 时，单次 Confirm 不得因 ErrOverloaded 提前置 FAILED；
// 容量恢复后 crash error_log 应成功入队并恰好上传一次。
func TestConfirmWaitsForCrashLogBeyondQueueWait(t *testing.T) {
	h := newFinishHarness(t, finishOptions{
		skipStore:   true,
		uploadQueue: 2,
		queueWait:   20 * time.Millisecond,
		blockUpload: true,
	})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("crash/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	select {
	case <-h.uploadEntered:
	case <-time.After(3 * time.Second):
		t.Fatal("uploader did not enter the sender")
	}

	saturated := false
	for step := 1; step <= 16; step++ {
		err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
			scalarRecord("crash/loss", int64(step), columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		})
		if err == nil {
			continue
		}
		if !errors.Is(err, ErrOverloaded) {
			t.Fatalf("fill upload queue: %v", err)
		}
		saturated = true
		break
	}
	if !saturated {
		t.Fatal("upload queue was never saturated")
	}

	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_CRASHED, "boom")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}

	type confirmOut struct {
		res ConfirmResult
		err error
	}
	done := make(chan confirmOut, 1)
	go func() {
		res, err := h.mgr.Confirm(ctx, h.id)
		done <- confirmOut{res: res, err: err}
	}()

	// 阻塞 80ms（4 倍于 queueWait=20ms），Confirm 必须仍在等待而非因 ErrOverloaded 失败退出。
	select {
	case out := <-done:
		t.Fatalf("Confirm returned prematurely while upload queue stayed full: %+v", out)
	case <-time.After(80 * time.Millisecond):
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining {
		t.Fatalf("state while waiting = %s, want draining", state)
	}

	close(h.gate)

	select {
	case out := <-done:
		if out.err != nil || out.res.ReleaseFailed {
			t.Fatalf("Confirm = %+v, want success", out)
		}
	case <-time.After(3 * time.Second):
		t.Fatal("Confirm did not complete after upload capacity freed")
	}

	if h.session.hasPendingLog() {
		t.Fatal("crash error log must be cleared after upload")
	}
	var logs []*recordv1.Record
	for _, rec := range h.uploads.snapshot() {
		if rec.GetLog() != nil {
			logs = append(logs, rec)
		}
	}
	if len(logs) != 1 || logs[0].GetLog().GetLine() != "boom" {
		t.Fatalf("uploaded logs = %v, want single 'boom' log", logs)
	}
	if _, err := h.mgr.State(h.id); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

func TestDeliverRunFinishConcurrentCallsWriteOnce(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	const callers = 4
	var wg sync.WaitGroup
	errs := make(chan error, callers)
	for range callers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			errs <- h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
		}()
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		if err != nil {
			t.Fatalf("concurrent DeliverRunFinish: %v", err)
		}
	}
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	waitDrained(t, h.session)

	finishCount := 0
	for _, item := range readStore(t, h.runFile) {
		if item.GetNum() == record.NumFinish {
			finishCount++
		}
	}
	if finishCount != 1 {
		t.Fatalf("stored %d finish records, want 1", finishCount)
	}
}

func TestDeliverRunFinishRejectsNilRecord(t *testing.T) {
	m := newTestManager()
	handle := installBare(t, m, sessionOptions{})
	if err := m.DeliverRunFinish(context.Background(), handle, nil); !errors.Is(err, ErrInvalidFinishRecord) {
		t.Fatalf("nil finish record err=%v, want ErrInvalidFinishRecord", err)
	}
	if err := m.DeliverRunFinish(context.Background(), "", finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); !errors.Is(err, ErrEmptySessionID) {
		t.Fatalf("empty session err=%v, want ErrEmptySessionID", err)
	}
}

func TestDeliverRunFinishSkipStore(t *testing.T) {
	h := newFinishHarness(t, finishOptions{skipStore: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertSaves(ctx, h.id, []*savev1.SaveRecord{{
		Name:       "config.yaml",
		Policy:     savev1.SavePolicy_SAVE_POLICY_NOW,
		Type:       savev1.SaveType_SAVE_TYPE_CONFIG,
		Payload:    []byte("key: value\n"),
		SourcePath: "",
	}}); err != nil {
		t.Fatalf("UpsertSaves: %v", err)
	}
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining && state != RunStateDrained {
		t.Fatalf("state after DeliverRunFinish = %s, want draining or drained", state)
	}
	// skip_store 下 finishCommitted 同样置位，Confirm 完成通道关闭、排空与释放。
	if _, err := h.mgr.Confirm(ctx, h.id); err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if got := len(h.uploads.snapshot()); got != 1 {
		t.Fatalf("uploaded %d records, want only the config save", got)
	}
	if _, err := h.mgr.State(h.id); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
}

// Shutdown 抢占收尾尾部时不写 FinishRecord：此时 DeliverRunFinish 不能假成功，
// 否则会话会带着缺失的 -2 结束。
func TestDeliverRunFinishAfterShutdownReturnsFailure(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.session.beginDrain(); err != nil {
		t.Fatalf("beginDrain: %v", err)
	}
	err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
	if !errors.Is(err, ErrRunFailed) {
		t.Fatalf("DeliverRunFinish err=%v, want ErrRunFailed", err)
	}
	if h.session.finishCommitted {
		t.Fatal("shutdown must not mark the finish record as committed")
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining {
		t.Fatalf("state = %s, want draining", state)
	}
}

// 对照：正常收尾排空后重复投递仍按幂等成功返回。
func TestDeliverRunFinishRepeatsAfterDrain(t *testing.T) {
	h := newFinishHarness(t, finishOptions{})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	waitDrained(t, h.session)
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("repeated DeliverRunFinish after drain: %v", err)
	}
}

func TestRunStateFromProto(t *testing.T) {
	cases := []struct {
		in   runv1.RunState
		want api.RunState
	}{
		// RUNNING 不是合法的收尾状态，Python 侧同样按正常结束上报。
		{runv1.RunState_RUN_STATE_RUNNING, api.RunStateFinished},
		{runv1.RunState_RUN_STATE_FINISHED, api.RunStateFinished},
		{runv1.RunState_RUN_STATE_CRASHED, api.RunStateCrashed},
		{runv1.RunState_RUN_STATE_ABORTED, api.RunStateAborted},
		{runv1.RunState(99), api.RunStateFinished},
	}
	for _, tc := range cases {
		if got := runStateFromProto(tc.in); got != tc.want {
			t.Fatalf("runStateFromProto(%v) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

func newCheckpointTestSession(
	t *testing.T,
	mode config.Mode,
	sender uploader.Sender,
	cpInterval time.Duration,
	onTick func(),
) (*Session, string, string) {
	t.Helper()
	cfg := config.Defaults()
	cfg.Mode = mode
	cfg.RunID = "cp-run"
	cfg.RunDir = t.TempDir()
	cfg.RecordBatch = 16
	cfg.RecordInterval = 0.005

	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	runFile, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}
	cpFile, err := run.CheckpointFile()
	if err != nil {
		t.Fatalf("CheckpointFile: %v", err)
	}
	writer, err := store.Create(runFile)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{
		Queue:         64,
		FlushBytes:    1,
		FlushInterval: 2 * time.Millisecond,
		IndexStride:   1,
	})
	if err != nil {
		_ = writer.Close()
		t.Fatalf("NewAsyncWriter: %v", err)
	}

	b := budget.New(budget.Limits{Records: 64, Bytes: 1 << 20})
	tracker := uploader.NewTracker()
	var trackerProgress *progress.Tracker
	var onStored func(*recordv1.Record, uint64)
	var watermark func() uint64
	if mode == config.ModeOnline {
		trackerProgress = progress.NewTracker(0)
		onStored = trackerProgress.OnStored
		watermark = trackerProgress.Watermark
	}

	intake, err := buffer.New(buffer.Options{
		StoreOnly:   mode != config.ModeOnline,
		Store:       async.In(),
		UploadQueue: 64,
		Budget:      b,
		Flushed: func() uint64 {
			seq, _ := async.Flushed()
			return seq
		},
		OnAdmit:  tracker.AddTotalRecords,
		OnStored: onStored,
	})
	if err != nil {
		_ = async.Close()
		t.Fatalf("buffer.New: %v", err)
	}

	s := newSession(sessionOptions{})
	s.cfg = cfg
	s.run = run
	s.store = async
	s.intake = intake
	s.budget = b
	s.tracker = tracker
	s.metrics = metrics.New()
	s.counter = counter.New(0)
	s.epoch = counter.New(0)
	s.watcher = watcher.NullWatcher{}
	s.runFile = runFile
	s.runFileCreated = true
	s.progress = trackerProgress
	if mode == config.ModeOnline && trackerProgress != nil {
		s.checkpoint = checkpoint.New(checkpoint.Options{
			File:      cpFile,
			RunID:     run.RunID,
			RunFile:   runFile,
			Store:     async,
			Watermark: watermark,
			Interval:  cpInterval,
			OnTick:    onTick,
		})
	}

	if mode == config.ModeOnline {
		up, err := uploader.New(uploader.Options{
			In:     intake.UploadOut(),
			Budget: b,
			Senders: map[record.Type]uploader.Sender{
				record.TypeScalar: sender,
			},
			Batch:     16,
			Interval:  5 * time.Millisecond,
			Tracker:   tracker,
			OnConfirm: trackerProgress.Resolve,
			OnDrop:    trackerProgress.ResolveOne,
			OnDrained: func() {
				s.checkpoint.Stop()
				s.checkpoint.Flush()
				s.markDrained()
			},
		})
		if err != nil {
			_ = async.Close()
			t.Fatalf("uploader.New: %v", err)
		}
		s.uploader = up
		s.drainFn = up.Wait
	}
	return s, runFile, cpFile
}

func waitStoreFlushed(t *testing.T, s *store.AsyncWriter, want uint64) {
	t.Helper()
	deadline := time.Now().Add(2 * time.Second)
	for {
		seq, _ := s.Flushed()
		if seq >= want {
			return
		}
		if time.Now().After(deadline) {
			t.Fatalf("flushed seq = %d, want >= %d", seq, want)
		}
		time.Sleep(time.Millisecond)
	}
}

func commitToIntake(t *testing.T, in *buffer.Intake, cfg config.Settings, recs ...*recordv1.Record) {
	t.Helper()
	tk, err := in.Reserve(context.Background(), planner.DemandFor(recs, cfg))
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := in.Commit(tk, recs); err != nil {
		t.Fatalf("Commit: %v", err)
	}
}

func TestCheckpointBoundsAndFinalFlushOnline(t *testing.T) {
	gate := make(chan struct{})
	entered := make(chan struct{})
	collect := &collectingSender{notify: make(chan struct{}, 8)}
	gated := &gatedSender{inner: collect, gate: gate, entered: entered}

	ticked := make(chan struct{}, 16)
	s, runFile, cpFile := newCheckpointTestSession(
		t, config.ModeOnline, gated, 5*time.Millisecond,
		func() {
			select {
			case ticked <- struct{}{}:
			default:
			}
		},
	)
	t.Cleanup(func() {
		s.stopOnlineDrain()
		s.checkpoint.Stop()
		if s.uploader != nil {
			_ = s.uploader.Close()
		}
		if s.store != nil {
			_ = s.store.Close()
		}
	})
	s.checkpoint.Start()

	// seq=1: Start (立即完结); seq=2: Scalar (将被 gatedSender 阻塞上传)。
	startRec := backfillStart()
	scRec := backfillScalar(1)
	commitToIntake(t, s.intake, s.cfg, startRec, scRec)

	waitStoreFlushed(t, s.store, 2)
	<-entered // uploader 已拿到 seq=2 并阻塞在 gate，此时 watermark=1 < flushedSeq=2。

	select {
	case <-ticked:
	case <-time.After(2 * time.Second):
		t.Fatal("timeout waiting for initial checkpoint tick")
	}

	cp, ok, err := store.LoadCheckpoint(cpFile, s.run.RunID)
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint(blocked) = (%v, %v), want (nil, true)", err, ok)
	}
	if cp.StoreSeq != 1 {
		t.Fatalf("StoreSeq while seq=2 unconfirmed = %d, want 1", cp.StoreSeq)
	}
	flushedSeq, flushedOff := s.store.Flushed()
	if flushedSeq != 2 {
		t.Fatalf("flushedSeq = %d, want 2", flushedSeq)
	}
	wantOff, err := store.EndOffset(runFile, s.store.Index(), 1, flushedOff)
	if err != nil {
		t.Fatalf("EndOffset: %v", err)
	}
	if cp.Offset != wantOff {
		t.Fatalf("Offset = %d, want %d", cp.Offset, wantOff)
	}

	// 放行上传并触发收尾：OnDrained 应停止协程并落盘最终连续水位（含 seq=3 Finish）。
	close(gate)
	collect.waitFor(t, 1)

	if err = s.finish(context.Background(), &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
		t.Fatalf("finish: %v", err)
	}
	if err = s.waitDrained(context.Background()); err != nil {
		t.Fatalf("waitDrained: %v", err)
	}

	cpFinal, ok, err := store.LoadCheckpoint(cpFile, s.run.RunID)
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint(final) = (%v, %v), want (nil, true)", err, ok)
	}
	finalSeq, finalOff := s.store.Flushed()
	if finalSeq != 3 || cpFinal.StoreSeq != finalSeq || cpFinal.Offset != finalOff {
		t.Fatalf("final cp = %+v, want seq=%d offset=%d", cpFinal, finalSeq, finalOff)
	}
}

func TestCheckpointLocalAndOfflineWriteNoFile(t *testing.T) {
	for _, mode := range []config.Mode{config.ModeLocal, config.ModeOffline} {
		t.Run(string(mode), func(t *testing.T) {
			s, _, cpFile := newCheckpointTestSession(t, mode, nil, 0, nil)
			t.Cleanup(func() {
				s.checkpoint.Stop()
				if s.store != nil {
					_ = s.store.Close()
				}
			})
			s.checkpoint.Start()

			commitToIntake(t, s.intake, s.cfg, backfillStart(), backfillScalar(1), backfillScalar(2))
			if err := s.finish(context.Background(), &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}); err != nil {
				t.Fatalf("finish: %v", err)
			}

			if _, err := os.Stat(cpFile); !os.IsNotExist(err) {
				t.Fatalf("%s mode must not write checkpoint file %s, got stat err=%v", mode, cpFile, err)
			}
		})
	}
}

func commitBackfill(t *testing.T, s *Session, records ...*recordv1.Record) {
	t.Helper()
	tk, err := s.intake.Reserve(context.Background(), planner.DemandFor(records, s.cfg))
	if err != nil {
		t.Fatalf("Reserve: %v", err)
	}
	if err := s.intake.Commit(tk, records); err != nil {
		t.Fatalf("Commit: %v", err)
	}
}

func backfillNums(records []*recordv1.Record) []int64 {
	nums := make([]int64, len(records))
	for i, rec := range records {
		nums[i] = rec.GetNum()
	}
	return nums
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
	if err := h.session.backfill.RunOnce(ctx); err != nil {
		t.Fatalf("RunOnce: %v", err)
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
	if err := h.session.backfill.RunOnce(ctx); err != nil {
		t.Fatalf("RunOnce: %v", err)
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

func TestFinishRejectsDegradedStoreDuringPausedBackfill(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 1})
	ctx, cancel := context.WithTimeout(context.Background(), time.Second)
	defer cancel()

	// 单批两条记录超过预算，首次提交即发布暂停。
	commitBackfill(t, h.session, backfillScalar(1), backfillScalar(2))
	if paused, seq := h.session.intake.PauseState(); !paused || seq != 1 {
		t.Fatalf("PauseState = (%v, %d), want (true, 1)", paused, seq)
	}
	h.session.storeDegraded.Store(true)

	// DeliverRunFinish 只负责本地收尾；暂停区间补投在 Confirm 时发现本地已损坏并拒绝上报终态。
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	_, err := h.mgr.Confirm(ctx, h.id)
	if !errors.Is(err, ErrRunFailed) {
		t.Fatalf("Confirm err = %v, want ErrRunFailed", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateFailed {
		t.Fatalf("state = %s, want failed", state)
	}
}

// 站点降级期间调用 DeliverRunFinish 必须快速完成本地收尾；Confirm 阻塞等待站点恢复，
// 恢复后自动从 run-*.swanlab 补投并排空上传，再执行 releaseFn 释放会话。
func TestDeliverRunFinishSucceedsWhileSiteSuspendedAndConfirmWaitsForRecovery(t *testing.T) {
	var houseHealthy atomic.Bool
	var released atomic.Bool
	reg := health.NewRegistry(nil)
	h := newFinishHarness(t, finishOptions{
		budgetRecords: 10,
		healthReg:     reg,
		healthOpts: &health.SiteOptions{
			Key:                  "outage-finish",
			Interval:             10 * time.Millisecond,
			Timeout:              50 * time.Millisecond,
			DegradeThreshold:     1,
			RecoverThreshold:     1,
			InitialHouseDegraded: true,
			Server:               func(context.Context) bool { return true },
			House:                func(context.Context) bool { return houseHealthy.Load() },
		},
		releaseFn: func(context.Context) error {
			released.Store(true)
			return nil
		},
	})

	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("outage/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("outage/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("outage/loss", 2, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	if paused, seq := h.session.intake.PauseState(); !paused || seq != 1 {
		t.Fatalf("PauseState = (%v, %d), want (true, 1)", paused, seq)
	}
	if got := len(h.uploads.snapshot()); got != 0 {
		t.Fatalf("uploaded %d records while site suspended, want 0", got)
	}

	// 站点仍然降级，但 DeliverRunFinish 只做本地收尾，200ms 短超时内必须立即成功。
	deliverCtx, deliverCancel := context.WithTimeout(ctx, 200*time.Millisecond)
	defer deliverCancel()
	if err := h.mgr.DeliverRunFinish(deliverCtx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish while site suspended: %v", err)
	}
	if state := h.session.stateSnapshot(); state != RunStateDraining {
		t.Fatalf("state after DeliverRunFinish = %s, want draining", state)
	}

	// Confirm 必须一直阻塞等待站点恢复，不能因站点降级提前失败退出。
	confirmDone := make(chan error, 1)
	go func() {
		_, err := h.mgr.Confirm(ctx, h.id)
		confirmDone <- err
	}()

	select {
	case err := <-confirmDone:
		t.Fatalf("Confirm returned prematurely while site still degraded: %v", err)
	case <-time.After(60 * time.Millisecond):
	}
	if released.Load() {
		t.Fatal("releaseFn must not run before backfill & drain complete")
	}

	// 恢复站点健康：backfill.Complete 唤醒补投、uploader 排空、stop_experiment 释放会话。
	houseHealthy.Store(true)
	select {
	case err := <-confirmDone:
		if err != nil {
			t.Fatalf("Confirm after site recovery: %v", err)
		}
	case <-time.After(3 * time.Second):
		t.Fatal("Confirm did not complete after site recovery")
	}

	if !released.Load() {
		t.Fatal("releaseFn was not called")
	}
	if _, err := h.mgr.State(h.id); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("State after confirm err=%v, want ErrSessionNotFound", err)
	}
	var scalars int
	for _, rec := range h.uploads.snapshot() {
		if record.TypeOf(rec) == record.TypeScalar {
			scalars++
		}
	}
	if scalars != 3 {
		t.Fatalf("uploaded scalars = %d, want 3", scalars)
	}
	if stored := readStore(t, h.runFile); len(stored) != 5 {
		t.Fatalf("store has %d records, want 1 column + 3 scalars + 1 finish", len(stored))
	}
}

// 补投协程由暂停通知唤醒，恢复后自行停止。
func TestBackfillSupervisorResumesAfterPause(t *testing.T) {
	h := newFinishHarness(t, finishOptions{budgetRecords: 2, blockUpload: true})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	h.session.backfill.Start()
	commitBackfill(t, h.session, backfillScalar(1))
	commitBackfill(t, h.session, backfillScalar(2))
	commitBackfill(t, h.session, backfillScalar(3))
	commitBackfill(t, h.session, backfillScalar(4))
	// 暂停期间的 END 保存不占预算，后台补投仍能推进。
	commitBackfill(t, h.session, backfillEndSave(5))
	close(h.gate)
	h.session.backfill.Wake()

	// 4 条标量送达说明后台补投已在推进（前两条来自实时通道）。
	h.uploads.waitFor(t, 4)
	// 交接完成前取消可能停在补投之后，这里显式 join 并补完剩余一步。
	h.session.backfill.Stop()
	if err := h.session.backfill.RunOnce(ctx); err != nil {
		t.Fatalf("complete backfill: %v", err)
	}
	// 暂停已清除：新提交的标量直接进入上传通道。
	commitBackfill(t, h.session, backfillScalar(6))
	h.uploads.waitFor(t, 5)

	if err := h.session.intake.Close(ctx); err != nil {
		t.Fatalf("intake.Close: %v", err)
	}
	h.session.uploader.RequestFinish()
	if err := h.session.uploader.Wait(ctx); err != nil {
		t.Fatalf("uploader.Wait: %v", err)
	}
	uploads := h.uploads.snapshot()
	if len(uploads) != 6 {
		t.Fatalf("uploaded %d records, want 6 (5 scalars + END save)", len(uploads))
	}
	if last := uploads[len(uploads)-1]; record.TypeOf(last) != record.TypeSave {
		t.Fatalf("last upload = %v, want the END save", record.TypeOf(last))
	}
	if paused, seq := h.session.intake.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState = (%v, %d), want resumed", paused, seq)
	}
}

// 仅调用 DeliverRunFinish、不提前调 Confirm 时，online 会话也必须自行推进到 DRAINED 与 CORE_STATE_FINISHED，
// 否则 Python 进度轮询会一直等待 CORE_STATE_FINISHED 而双方互锁。
func TestDeliverRunFinishDrainsOnlineWithoutConfirm(t *testing.T) {
	t.Run("healthy", func(t *testing.T) {
		h := newFinishHarness(t, finishOptions{})
		h.session.tracker.SetState(operationv1.CoreState_CORE_STATE_RUNNING)
		ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
		defer cancel()

		if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
			scalarRecord("drain/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
			scalarRecord("drain/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		}); err != nil {
			t.Fatalf("UpsertScalars: %v", err)
		}
		if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
			t.Fatalf("DeliverRunFinish: %v", err)
		}
		waitDrained(t, h.session)

		if state := h.session.stateSnapshot(); state != RunStateDrained {
			t.Fatalf("state = %s, want drained", state)
		}
		stats, err := h.mgr.Stats(h.id)
		if err != nil {
			t.Fatalf("Stats: %v", err)
		}
		if stats.GetState() != operationv1.CoreState_CORE_STATE_FINISHED {
			t.Fatalf("stats state = %v, want CORE_STATE_FINISHED", stats.GetState())
		}
		var scalars int
		for _, rec := range h.uploads.snapshot() {
			if record.TypeOf(rec) == record.TypeScalar {
				scalars++
			}
		}
		if scalars != 2 {
			t.Fatalf("uploaded scalars = %d, want 2", scalars)
		}
	})

	t.Run("suspended_then_recovered", func(t *testing.T) {
		var houseHealthy atomic.Bool
		reg := health.NewRegistry(nil)
		h := newFinishHarness(t, finishOptions{
			budgetRecords: 10,
			healthReg:     reg,
			healthOpts: &health.SiteOptions{
				Key:                  "no-confirm-outage",
				Interval:             10 * time.Millisecond,
				Timeout:              50 * time.Millisecond,
				DegradeThreshold:     1,
				RecoverThreshold:     1,
				InitialHouseDegraded: true,
				Server:               func(context.Context) bool { return true },
				House:                func(context.Context) bool { return houseHealthy.Load() },
			},
		})
		h.session.tracker.SetState(operationv1.CoreState_CORE_STATE_RUNNING)
		ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
		defer cancel()

		if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
			scalarRecord("drain/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
			scalarRecord("drain/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
			scalarRecord("drain/loss", 2, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		}); err != nil {
			t.Fatalf("UpsertScalars: %v", err)
		}
		if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
			t.Fatalf("DeliverRunFinish: %v", err)
		}

		// 站点恢复前不得提前标为 FINISHED。
		time.Sleep(40 * time.Millisecond)
		stats, err := h.mgr.Stats(h.id)
		if err != nil {
			t.Fatalf("Stats: %v", err)
		}
		if stats.GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
			t.Fatalf("stats state while suspended = %v, want CORE_STATE_RUNNING", stats.GetState())
		}

		houseHealthy.Store(true)
		waitDrained(t, h.session)

		if state := h.session.stateSnapshot(); state != RunStateDrained {
			t.Fatalf("state after recovery = %s, want drained", state)
		}
		stats, err = h.mgr.Stats(h.id)
		if err != nil {
			t.Fatalf("Stats after recovery: %v", err)
		}
		if stats.GetState() != operationv1.CoreState_CORE_STATE_FINISHED {
			t.Fatalf("stats state after recovery = %v, want CORE_STATE_FINISHED", stats.GetState())
		}
		var scalars int
		for _, rec := range h.uploads.snapshot() {
			if record.TypeOf(rec) == record.TypeScalar {
				scalars++
			}
		}
		if scalars != 3 {
			t.Fatalf("uploaded scalars = %d, want 3", scalars)
		}
	})
}

// 后台在线排空长期阻塞在 finishMu 时，Shutdown 必须能立即打断它而不是等满整段退避/退出窗口。
func TestShutdownInterruptsBackgroundOnlineDrain(t *testing.T) {
	reg := health.NewRegistry(nil)
	h := newFinishHarness(t, finishOptions{
		budgetRecords: 10,
		healthReg:     reg,
		healthOpts: &health.SiteOptions{
			Key:                  "shutdown-interrupt",
			Interval:             10 * time.Millisecond,
			Timeout:              50 * time.Millisecond,
			DegradeThreshold:     1,
			RecoverThreshold:     1,
			InitialHouseDegraded: true,
			Server:               func(context.Context) bool { return true },
			House:                func(context.Context) bool { return false },
		},
	})
	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()

	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("shutdown/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("shutdown/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	if err := h.mgr.DeliverRunFinish(ctx, h.id, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	time.Sleep(30 * time.Millisecond)

	start := time.Now()
	shutdownCtx, shutdownCancel := context.WithTimeout(ctx, 2*time.Second)
	defer shutdownCancel()
	if err := h.mgr.Shutdown(shutdownCtx); err != nil {
		t.Fatalf("Shutdown: %v", err)
	}
	if elapsed := time.Since(start); elapsed >= 500*time.Millisecond {
		t.Fatalf("Shutdown took %v, expected immediate interruption of background drain", elapsed)
	}
}

// 上传通道持续满时，shutdown 对崩溃 error_log 的补投必须受短预算限制，不能把整段 ctx 耗光后再去关闭通道。
func TestShutdownBoundsCrashLogWait(t *testing.T) {
	origBudget := shutdownCrashLogBudget
	shutdownCrashLogBudget = 80 * time.Millisecond
	t.Cleanup(func() { shutdownCrashLogBudget = origBudget })

	cfg := config.Defaults()
	cfg.Mode = config.ModeOnline
	cfg.SkipStore = true
	cfg.RunID = "shutdown-crash-log"

	b := budget.New(budget.Limits{Records: 8, Bytes: 1 << 20})
	intake, err := buffer.New(buffer.Options{
		SkipStore:   true,
		UploadQueue: 1,
		Budget:      b,
		QueueWait:   20 * time.Millisecond,
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}

	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	t.Cleanup(func() {
		s.stopOnlineDrain()
		s.cancel()
	})
	s.cfg = cfg
	s.intake = intake
	s.budget = b
	s.counter = counter.New(0)
	s.epoch = counter.New(0)
	s.prepareFinish(finishRecord(runv1.RunState_RUN_STATE_CRASHED, "boom"))

	// 占满唯一的上传通道槽位且不启动 uploader，确保 commitCrashLog 持续遇到 ErrOverloaded。
	if err = intake.CommitUploadOnly(context.Background(), backfillScalar(1)); err != nil {
		t.Fatalf("fill upload queue: %v", err)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 3*time.Second)
	defer cancel()
	start := time.Now()
	err = s.shutdown(ctx)
	elapsed := time.Since(start)
	if !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("shutdown err = %v, want DeadlineExceeded from bounded crash log budget", err)
	}
	if elapsed >= 500*time.Millisecond {
		t.Fatalf("shutdown took %v, exceeded bounded crash log budget", elapsed)
	}
	if ctx.Err() != nil {
		t.Fatalf("outer shutdown ctx expired prematurely: %v", ctx.Err())
	}
}
