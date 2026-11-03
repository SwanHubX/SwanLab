package manager

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

type upsertOptions struct {
	skipStore bool
	batch     int32
	queue     int
	section   int32
}

type upsertHarness struct {
	mgr     *Manager
	id      string
	session *Session
	storeCh chan store.Request
	upload  <-chan *recordv1.Record
}

func newUpsertHarness(t *testing.T, opts upsertOptions) *upsertHarness {
	t.Helper()
	if opts.batch <= 0 {
		opts.batch = 16
	}
	if opts.queue <= 0 {
		opts.queue = queueSize(opts.batch)
	}
	cfg := config.Defaults()
	cfg.RunID = "upsert"
	cfg.RunDir = t.TempDir()
	cfg.RecordBatch = opts.batch
	cfg.RecordInterval = 1
	cfg.SectionRule = opts.section
	cfg.SkipStore = opts.skipStore

	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.run = run
	s.metrics = metrics.New()
	s.counter = counter.New(0)
	s.epoch = counter.New(0)
	s.watcher = watcher.NullWatcher{}
	s.tracker = uploader.NewTracker()

	var storeCh chan store.Request
	if !opts.skipStore {
		storeCh = make(chan store.Request, opts.queue)
	}
	intake, err := buffer.New(buffer.Options{
		SkipStore:   opts.skipStore,
		Store:       storeCh,
		UploadQueue: opts.queue,
		Budget:      budget.New(budget.Limits{}),
		QueueWait:   20 * time.Millisecond,
		Flushed:     func() uint64 { return 0 },
		OnAdmit:     s.tracker.AddTotalRecords,
	})
	if err != nil {
		t.Fatalf("buffer.New: %v", err)
	}
	s.intake = intake

	mgr := newTestManager()
	id, err := mgr.install(s)
	if err != nil {
		t.Fatalf("install: %v", err)
	}
	t.Cleanup(func() {
		ctx, cancel := context.WithTimeout(context.Background(), time.Second)
		defer cancel()
		_ = intake.Close(ctx)
		if s.cancel != nil {
			s.cancel()
		}
	})
	return &upsertHarness{mgr: mgr, id: id, session: s, storeCh: storeCh, upload: intake.UploadOut()}
}

func drainRecords(ch <-chan *recordv1.Record) []*recordv1.Record {
	var out []*recordv1.Record
	for {
		select {
		case rec, ok := <-ch:
			if !ok {
				return out
			}
			out = append(out, rec)
		default:
			return out
		}
	}
}

func drainStore(ch chan store.Request) []store.Request {
	var out []store.Request
	for {
		select {
		case req := <-ch:
			out = append(out, req)
		default:
			return out
		}
	}
}

func scalarRecord(key string, step int64, typ columnv1.ColumnType) *datav1.ScalarRecord {
	return &datav1.ScalarRecord{Key: key, Step: step, Type: typ}
}

func TestUpsertColumnsSkipsDefinedAndDuplicates(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()
	columns := []*columnv1.ColumnRecord{
		{ColumnKey: "train/loss", ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR},
		{ColumnKey: "train/loss", ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR},
		{ColumnKey: "train/image", ColumnType: columnv1.ColumnType_COLUMN_TYPE_IMAGE},
	}
	if err := h.mgr.UpsertColumns(ctx, h.id, columns); err != nil {
		t.Fatalf("UpsertColumns: %v", err)
	}
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 2 {
		t.Fatalf("uploaded %d column records, want 2", len(uploaded))
	}
	if got := uploaded[0].GetColumn().GetColumnKey(); got != "train/loss" {
		t.Fatalf("first column = %q, want train/loss", got)
	}
	if got := uploaded[1].GetColumn().GetColumnKey(); got != "train/image" {
		t.Fatalf("second column = %q, want train/image", got)
	}
	if !h.session.metrics.Defined("train/loss") || !h.session.metrics.Defined("train/image") {
		t.Fatal("columns were not defined in metrics")
	}
	if minStep, ok := h.session.metrics.MinStep("train/loss"); !ok || minStep != -1 {
		t.Fatalf("min step = %d (defined=%v), want -1", minStep, ok)
	}

	if err := h.mgr.UpsertColumns(ctx, h.id, columns); err != nil {
		t.Fatalf("repeated UpsertColumns: %v", err)
	}
	if got := drainRecords(h.upload); len(got) != 0 {
		t.Fatalf("repeated columns uploaded %d records, want 0", len(got))
	}
	if got := len(drainStore(h.storeCh)); got != 2 {
		t.Fatalf("store received %d records, want 2", got)
	}
}

func TestUpsertScalarsAutoColumnAndStepDedup(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()
	records := []*datav1.ScalarRecord{
		scalarRecord("train/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("train/loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("train/loss", 1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("train/loss", -1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}
	if err := h.mgr.UpsertScalars(ctx, h.id, records); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 3 {
		t.Fatalf("uploaded %d records, want column + 2 scalars", len(uploaded))
	}
	column := uploaded[0].GetColumn()
	if column.GetColumnKey() != "train/loss" || column.GetSectionName() != "train" {
		t.Fatalf("auto column = %+v, want key train/loss in section train", column)
	}
	if got := uploaded[1].GetScalar().GetStep(); got != 0 {
		t.Fatalf("first scalar step = %d, want 0", got)
	}
	if got := uploaded[2].GetScalar().GetStep(); got != 1 {
		t.Fatalf("second scalar step = %d, want 1", got)
	}
	stats, err := h.mgr.Stats(h.id)
	if err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if stats.GetTotalRecords() != 3 {
		t.Fatalf("total records = %d, want 3 (column + two scalars)", stats.GetTotalRecords())
	}

	if err := h.mgr.UpsertScalars(ctx, h.id, records[:1]); err != nil {
		t.Fatalf("duplicate UpsertScalars: %v", err)
	}
	if got := drainRecords(h.upload); len(got) != 0 {
		t.Fatalf("duplicate step uploaded %d records, want 0", len(got))
	}
	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("train/loss", 2, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("next step UpsertScalars: %v", err)
	}
	if got := drainRecords(h.upload); len(got) != 1 {
		t.Fatalf("next step uploaded %d records, want 1", len(got))
	}
}

func TestUpsertScalarsSkipsDeterministicErrors(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()
	columns := []*columnv1.ColumnRecord{{ColumnKey: "loss", ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR}}
	if err := h.mgr.UpsertColumns(ctx, h.id, columns); err != nil {
		t.Fatalf("UpsertColumns: %v", err)
	}
	drainRecords(h.upload)

	records := []*datav1.ScalarRecord{
		scalarRecord("loss", metrics.MaxStep+1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("loss", metrics.MaxStep, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("loss", 0, columnv1.ColumnType_COLUMN_TYPE_IMAGE),
	}
	if err := h.mgr.UpsertScalars(ctx, h.id, records); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 1 {
		t.Fatalf("uploaded %d records, want only the max-step record", len(uploaded))
	}
	if got := uploaded[0].GetScalar().GetStep(); got != metrics.MaxStep {
		t.Fatalf("accepted step = %d, want %d", got, metrics.MaxStep)
	}

	// MaxStep 已被接纳：重复提交必须跳过；超限值再次提交也不能产出记录。
	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{
		scalarRecord("loss", metrics.MaxStep, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("loss", metrics.MaxStep+1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}); err != nil {
		t.Fatalf("boundary UpsertScalars: %v", err)
	}
	if got := drainRecords(h.upload); len(got) != 0 {
		t.Fatalf("boundary reprocess uploaded %d records, want 0", len(got))
	}
}

func TestUpsertScalarsKeepsBatchAfterRejectedRecord(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()
	records := []*datav1.ScalarRecord{
		scalarRecord("broken", metrics.MaxStep+1, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
		scalarRecord("ok", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR),
	}
	if err := h.mgr.UpsertScalars(ctx, h.id, records); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 3 {
		t.Fatalf("uploaded %d records, want two columns + ok scalar", len(uploaded))
	}
	if !h.session.metrics.Defined("broken") {
		t.Fatal("column should be defined before the step is rejected (Python parity)")
	}
	if !h.session.metrics.Defined("ok") {
		t.Fatal("valid record in the same batch must still be accepted")
	}
}

func TestUpsertMediaAutoColumnAndTypeMismatch(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()
	records := []*datav1.MediaRecord{
		{Key: "media/text", Step: 0, Type: columnv1.ColumnType_COLUMN_TYPE_TEXT},
		{Key: "media/text", Step: 0, Type: columnv1.ColumnType_COLUMN_TYPE_TEXT},
		{Key: "media/text", Step: 1, Type: columnv1.ColumnType_COLUMN_TYPE_IMAGE},
	}
	if err := h.mgr.UpsertMedia(ctx, h.id, records); err != nil {
		t.Fatalf("UpsertMedia: %v", err)
	}
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 2 {
		t.Fatalf("uploaded %d records, want column + one media record", len(uploaded))
	}
	column := uploaded[0].GetColumn()
	if column.GetColumnType() != columnv1.ColumnType_COLUMN_TYPE_TEXT {
		t.Fatalf("column type = %s, want TEXT", column.GetColumnType())
	}
	if len(column.GetColumnName()) != 0 {
		t.Fatalf("auto column should not invent a display name: %q", column.GetColumnName())
	}
	if got := uploaded[1].GetMedia().GetStep(); got != 0 {
		t.Fatalf("media step = %d, want 0", got)
	}
}

func TestUpsertLogsUsesResumeEpoch(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	h.session.epoch = counter.New(41)
	ctx := context.Background()
	logs := []*terminalv1.LogRecord{{Line: "first"}, {Line: "second"}}
	if err := h.mgr.UpsertLogs(ctx, h.id, logs); err != nil {
		t.Fatalf("UpsertLogs: %v", err)
	}
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 2 {
		t.Fatalf("uploaded %d logs, want 2", len(uploaded))
	}
	if uploaded[0].GetNum() != 1 || uploaded[1].GetNum() != 2 {
		t.Fatalf("log nums = %d,%d, want 1,2", uploaded[0].GetNum(), uploaded[1].GetNum())
	}
	if got := uploaded[0].GetLog().GetEpoch(); got != 42 {
		t.Fatalf("first epoch = %d, want 42", got)
	}
	if got := uploaded[1].GetLog().GetEpoch(); got != 43 {
		t.Fatalf("second epoch = %d, want 43", got)
	}
}

func TestUpsertSavesLinksAndDefersCustomEnd(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()
	sourceDir := t.TempDir()
	nowPath := filepath.Join(sourceDir, "now.txt")
	endPath := filepath.Join(sourceDir, "end.txt")
	for _, path := range []string{nowPath, endPath} {
		if err := os.WriteFile(path, []byte("data"), 0o644); err != nil {
			t.Fatalf("write %s: %v", path, err)
		}
	}
	nowSave := &savev1.SaveRecord{
		Name:       "now.txt",
		SourcePath: nowPath,
		Policy:     savev1.SavePolicy_SAVE_POLICY_NOW,
		Type:       savev1.SaveType_SAVE_TYPE_CUSTOM,
	}
	endSave := &savev1.SaveRecord{
		Name:       "end.txt",
		SourcePath: endPath,
		Policy:     savev1.SavePolicy_SAVE_POLICY_END,
		Type:       savev1.SaveType_SAVE_TYPE_CUSTOM,
	}
	configSave := &savev1.SaveRecord{
		Name:   "config.yaml",
		Policy: savev1.SavePolicy_SAVE_POLICY_NOW,
		Type:   savev1.SaveType_SAVE_TYPE_CONFIG,
	}
	saves := []*savev1.SaveRecord{nowSave, endSave, configSave}
	if err := h.mgr.UpsertSaves(ctx, h.id, saves); err != nil {
		t.Fatalf("UpsertSaves: %v", err)
	}

	filesDir := h.session.run.FilesDir()
	for _, name := range []string{"now.txt", "end.txt"} {
		if _, err := os.Stat(filepath.Join(filesDir, name)); err != nil {
			t.Fatalf("missing mirror link for %s: %v", name, err)
		}
	}
	if nowSave.GetTargetPath() == "" || endSave.GetTargetPath() == "" {
		t.Fatalf("target paths not filled: now=%q end=%q", nowSave.GetTargetPath(), endSave.GetTargetPath())
	}

	uploaded := drainRecords(h.upload)
	if len(uploaded) != 2 {
		t.Fatalf("uploaded %d saves before close, want internal + NOW", len(uploaded))
	}
	if got := uploaded[0].GetSave().GetType(); got != savev1.SaveType_SAVE_TYPE_CONFIG {
		t.Fatalf("first upload type = %s, want CONFIG", got)
	}
	if uploaded[0].GetNum() != -3 {
		t.Fatalf("config num = %d, want -3", uploaded[0].GetNum())
	}
	if got := uploaded[1].GetSave().GetPolicy(); got != savev1.SavePolicy_SAVE_POLICY_NOW {
		t.Fatalf("second upload policy = %s, want NOW", got)
	}
	if got := len(drainStore(h.storeCh)); got != 3 {
		t.Fatalf("store received %d records, want 3", got)
	}

	closeCtx, cancel := context.WithTimeout(context.Background(), time.Second)
	defer cancel()
	if err := h.session.intake.Close(closeCtx); err != nil {
		t.Fatalf("close intake: %v", err)
	}
	deferred := drainRecords(h.upload)
	if len(deferred) != 1 {
		t.Fatalf("deferred uploads = %d, want 1", len(deferred))
	}
	if got := deferred[0].GetSave().GetPolicy(); got != savev1.SavePolicy_SAVE_POLICY_END {
		t.Fatalf("deferred policy = %s, want END", got)
	}
	if deferred[0].GetNum() != 2 {
		t.Fatalf("deferred num = %d, want 2", deferred[0].GetNum())
	}
}

func TestUpsertSavesSkipStoreHasNoLinks(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{skipStore: true})
	ctx := context.Background()
	sourcePath := filepath.Join(t.TempDir(), "now.txt")
	if err := os.WriteFile(sourcePath, []byte("data"), 0o644); err != nil {
		t.Fatalf("write source: %v", err)
	}
	saves := []*savev1.SaveRecord{
		{
			Name:       "now.txt",
			SourcePath: sourcePath,
			Policy:     savev1.SavePolicy_SAVE_POLICY_NOW,
			Type:       savev1.SaveType_SAVE_TYPE_CUSTOM,
		},
		{
			Name:       "end.txt",
			SourcePath: sourcePath,
			Policy:     savev1.SavePolicy_SAVE_POLICY_END,
			Type:       savev1.SaveType_SAVE_TYPE_CUSTOM,
		},
	}
	if err := h.mgr.UpsertSaves(ctx, h.id, saves); err != nil {
		t.Fatalf("UpsertSaves: %v", err)
	}
	if _, err := os.Stat(h.session.run.FilesDir()); !os.IsNotExist(err) {
		t.Fatalf("skip_store must not create files dir, stat err=%v", err)
	}
	if saves[0].GetTargetPath() != "" || saves[1].GetTargetPath() != "" {
		t.Fatalf("skip_store must not fill target paths: %q %q", saves[0].GetTargetPath(), saves[1].GetTargetPath())
	}
	if got := drainRecords(h.upload); len(got) != 1 {
		t.Fatalf("uploaded %d saves before close, want only NOW", len(got))
	}
	closeCtx, cancel := context.WithTimeout(context.Background(), time.Second)
	defer cancel()
	if err := h.session.intake.Close(closeCtx); err != nil {
		t.Fatalf("close intake: %v", err)
	}
	if got := drainRecords(h.upload); len(got) != 1 {
		t.Fatalf("deferred uploads = %d, want 1", len(got))
	}
}

func TestUpsertRejectsWhenNotRunning(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()
	if err := drainBare(h.mgr, h.id); err != nil {
		t.Fatalf("Finish: %v", err)
	}
	calls := []struct {
		name string
		call func() error
	}{
		{"columns", func() error {
			return h.mgr.UpsertColumns(ctx, h.id, []*columnv1.ColumnRecord{{ColumnKey: "loss"}})
		}},
		{"scalars", func() error {
			return h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{scalarRecord("loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR)})
		}},
		{"media", func() error {
			return h.mgr.UpsertMedia(ctx, h.id, []*datav1.MediaRecord{{Key: "image", Type: columnv1.ColumnType_COLUMN_TYPE_IMAGE}})
		}},
		{"logs", func() error {
			return h.mgr.UpsertLogs(ctx, h.id, []*terminalv1.LogRecord{{Line: "line"}})
		}},
		{"saves", func() error {
			return h.mgr.UpsertSaves(ctx, h.id, []*savev1.SaveRecord{{Name: "file", Type: savev1.SaveType_SAVE_TYPE_CUSTOM}})
		}},
	}
	for _, tc := range calls {
		if err := tc.call(); !errors.Is(err, ErrNotRunning) {
			t.Fatalf("%s after drain err=%v, want ErrNotRunning", tc.name, err)
		}
	}
}

func TestUpsertReserveFailureLeavesMetricsUntouched(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{batch: 1, queue: 2})
	ctx := context.Background()
	col := &columnv1.ColumnRecord{ColumnKey: "loss", ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR}
	if err := h.session.metrics.DefineScalar("loss", col, -1); err != nil {
		t.Fatalf("DefineScalar: %v", err)
	}
	// 填满本地通道，强制 Reserve 在队列等待上限后失败。
	h.storeCh <- store.Request{}
	h.storeCh <- store.Request{}

	rec := scalarRecord("loss", 0, columnv1.ColumnType_COLUMN_TYPE_SCALAR)
	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{rec}); !errors.Is(err, ErrOverloaded) {
		t.Fatalf("UpsertScalars err=%v, want ErrOverloaded", err)
	}
	if got := drainRecords(h.upload); len(got) != 0 {
		t.Fatalf("failed reserve uploaded %d records, want 0", len(got))
	}

	// 位图未被污染：释放一个槽位后同一条 step 必须被重新接纳。
	<-h.storeCh
	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{rec}); err != nil {
		t.Fatalf("retry UpsertScalars: %v", err)
	}
	if got := drainRecords(h.upload); len(got) != 1 {
		t.Fatalf("retry uploaded %d records, want 1", len(got))
	}
}

func TestUpsertChunksByRecordBatch(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{batch: 2, queue: 32})
	ctx := context.Background()
	records := make([]*datav1.ScalarRecord, 0, 5)
	for step := int64(0); step < 5; step++ {
		records = append(records, scalarRecord("train/loss", step, columnv1.ColumnType_COLUMN_TYPE_SCALAR))
	}
	if err := h.mgr.UpsertScalars(ctx, h.id, records); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 6 {
		t.Fatalf("uploaded %d records, want column + 5 scalars", len(uploaded))
	}
	for i := 1; i < len(uploaded); i++ {
		if got := uploaded[i].GetScalar().GetStep(); got != int64(i-1) {
			t.Fatalf("record %d step = %d, want %d", i, got, i-1)
		}
	}
}

func TestUpsertRejectsInvalidRecords(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	ctx := context.Background()

	if err := h.mgr.UpsertColumns(ctx, h.id, []*columnv1.ColumnRecord{nil}); !errors.Is(err, ErrInvalidRecord) {
		t.Fatalf("nil column err=%v, want ErrInvalidRecord", err)
	}
	if err := h.mgr.UpsertColumns(ctx, h.id, []*columnv1.ColumnRecord{{}}); !errors.Is(err, ErrInvalidRecord) {
		t.Fatalf("empty column key err=%v, want ErrInvalidRecord", err)
	}
	if err := h.mgr.UpsertScalars(ctx, h.id, []*datav1.ScalarRecord{nil}); !errors.Is(err, ErrInvalidRecord) {
		t.Fatalf("nil scalar err=%v, want ErrInvalidRecord", err)
	}
	if err := h.mgr.UpsertMedia(ctx, h.id, []*datav1.MediaRecord{{}}); !errors.Is(err, ErrInvalidRecord) {
		t.Fatalf("empty media key err=%v, want ErrInvalidRecord", err)
	}
	if err := h.mgr.UpsertLogs(ctx, h.id, []*terminalv1.LogRecord{nil}); !errors.Is(err, ErrInvalidRecord) {
		t.Fatalf("nil log err=%v, want ErrInvalidRecord", err)
	}
	if err := h.mgr.UpsertSaves(ctx, h.id, []*savev1.SaveRecord{nil}); !errors.Is(err, ErrInvalidRecord) {
		t.Fatalf("nil save err=%v, want ErrInvalidRecord", err)
	}
	if err := h.mgr.UpsertSaves(ctx, h.id, []*savev1.SaveRecord{{Type: savev1.SaveType(99)}}); !errors.Is(err, ErrInvalidRecord) {
		t.Fatalf("unknown save type err=%v, want ErrInvalidRecord", err)
	}
	if err := h.mgr.UpsertScalars(ctx, "", nil); !errors.Is(err, ErrEmptySessionID) {
		t.Fatalf("empty session err=%v, want ErrEmptySessionID", err)
	}
	if err := h.mgr.UpsertScalars(ctx, "missing", nil); !errors.Is(err, ErrSessionNotFound) {
		t.Fatalf("unknown session err=%v, want ErrSessionNotFound", err)
	}
	// 空列表仍是 no-op。
	if err := h.mgr.UpsertColumns(ctx, h.id, []*columnv1.ColumnRecord{}); err != nil {
		t.Fatalf("empty columns: %v", err)
	}
}

func TestOnFileChangedReenqueuesLiveSave(t *testing.T) {
	h := newUpsertHarness(t, upsertOptions{})
	rec := &savev1.SaveRecord{
		Name:   "live.txt",
		Policy: savev1.SavePolicy_SAVE_POLICY_LIVE,
		Type:   savev1.SaveType_SAVE_TYPE_CUSTOM,
	}
	h.session.onFileChanged(rec)
	uploaded := drainRecords(h.upload)
	if len(uploaded) != 1 {
		t.Fatalf("uploaded %d records, want 1", len(uploaded))
	}
	if got := uploaded[0].GetSave().GetPolicy(); got != savev1.SavePolicy_SAVE_POLICY_LIVE {
		t.Fatalf("policy = %s, want LIVE", got)
	}

	if err := h.session.beginDrain(); err != nil {
		t.Fatalf("beginDrain: %v", err)
	}
	h.session.onFileChanged(rec)
	if got := drainRecords(h.upload); len(got) != 0 {
		t.Fatalf("callback after drain uploaded %d records, want 0", len(got))
	}
}
