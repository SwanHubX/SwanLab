// 本文件承载各类数据记录的预检与构建。
package manager

import (
	"context"
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/manager/planner"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func (s *Session) upsertColumns(ctx context.Context, columns []*columnv1.ColumnRecord) error {
	return commitChunks(ctx, s, columns, s.prepareColumns)
}

type columnDef struct {
	key    string
	column *columnv1.ColumnRecord
	kind   metrics.Kind
}

func (s *Session) prepareColumns(columns []*columnv1.ColumnRecord) (preparedBatch, error) {
	records := make([]*recordv1.Record, 0, len(columns))
	defs := make([]columnDef, 0, len(columns))
	seen := make(map[string]struct{}, len(columns))
	for _, column := range columns {
		if column == nil || column.GetColumnKey() == "" {
			return preparedBatch{}, fmt.Errorf("%w: column key is required", ErrInvalidRecord)
		}
		key := column.GetColumnKey()
		if _, duplicate := seen[key]; duplicate {
			logDefinedColumn(column)
			continue
		}
		seen[key] = struct{}{}
		if s.metrics.Defined(key) {
			logDefinedColumn(column)
			continue
		}
		defs = append(defs, columnDef{key: key, column: column, kind: columnKind(column)})
		records = append(records, builder.BuildColumnRecord(s.counter, column))
	}
	apply := func() error {
		for _, def := range defs {
			if err := s.metrics.Define(def.key, def.column, -1, def.kind); err != nil {
				return err
			}
		}
		return nil
	}
	return preparedBatch{records: records, apply: apply}, nil
}

func logDefinedColumn(column *columnv1.ColumnRecord) {
	if column.GetColumnClass() == columnv1.ColumnClass_COLUMN_CLASS_SYSTEM {
		console.Debugf("Column %s already has been defined, skipping", column.GetColumnKey())
		return
	}
	console.Warningf("Column %s already has been defined, skipping", column.GetColumnKey())
}

// columnKind 把列类型映射为指标类型：标量与直方图独立定义，其余按媒体列处理。
func columnKind(column *columnv1.ColumnRecord) metrics.Kind {
	switch column.GetColumnType() {
	case columnv1.ColumnType_COLUMN_TYPE_SCALAR:
		return metrics.KindScalar
	case columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM:
		return metrics.KindHistogram
	default:
		return metrics.KindMedia
	}
}

type metricRecordOps[T any] struct {
	kindLabel string
	key       func(T) string
	step      func(T) int64
	precheck  func(*planner.MetricPlanner, T) metrics.Precheck
	build     func(*counter.Counter, T) *recordv1.Record
	commit    func(T) error
}

func prepareMetricBatch[T any](s *Session, records []T, ops metricRecordOps[T]) (preparedBatch, error) {
	p := planner.NewMetricPlanner(s.metrics, s.cfg.SectionRule)
	planned := make([]*recordv1.Record, 0, len(records))
	accepted := make([]T, 0, len(records))
	for _, rec := range records {
		key := ops.key(rec)
		if key == "" {
			return preparedBatch{}, fmt.Errorf("%w: %s key is required", ErrInvalidRecord, ops.kindLabel)
		}
		verdict := ops.precheck(p, rec)
		if verdict.NeedsColumn {
			planned = append(planned, builder.BuildColumnRecord(s.counter, p.Column(key)))
		}
		switch {
		case verdict.Err != nil:
			console.Warningf("metric %q step %d skipped: %v", key, ops.step(rec), verdict.Err)
		case verdict.Decision == metrics.DecisionAccept:
			planned = append(planned, ops.build(s.counter, rec))
			accepted = append(accepted, rec)
		default:
			console.Debugf("Metric '%s' at step %d was skipped because it is duplicate or too old.",
				key, ops.step(rec))
		}
	}
	apply := func() error {
		if err := p.DefineAll(); err != nil {
			return err
		}
		for _, rec := range accepted {
			if err := ops.commit(rec); err != nil {
				return err
			}
		}
		return nil
	}
	return preparedBatch{records: planned, apply: apply}, nil
}

func (s *Session) upsertScalars(ctx context.Context, records []*datav1.ScalarRecord) error {
	return commitChunks(ctx, s, records, s.prepareScalars)
}

func (s *Session) prepareScalars(records []*datav1.ScalarRecord) (preparedBatch, error) {
	return prepareMetricBatch(s, records, metricRecordOps[*datav1.ScalarRecord]{
		kindLabel: "scalar",
		key: func(rec *datav1.ScalarRecord) string {
			if rec == nil {
				return ""
			}
			return rec.GetKey()
		},
		step:     func(rec *datav1.ScalarRecord) int64 { return rec.GetStep() },
		precheck: (*planner.MetricPlanner).PrecheckScalar,
		build:    builder.BuildScalarRecord,
		commit: func(rec *datav1.ScalarRecord) error {
			_, err := s.metrics.CommitScalar(rec, nil)
			return err
		},
	})
}

func (s *Session) upsertMedia(ctx context.Context, records []*datav1.MediaRecord) error {
	return commitChunks(ctx, s, records, s.prepareMedia)
}

func (s *Session) prepareMedia(records []*datav1.MediaRecord) (preparedBatch, error) {
	return prepareMetricBatch(s, records, metricRecordOps[*datav1.MediaRecord]{
		kindLabel: "media",
		key: func(rec *datav1.MediaRecord) string {
			if rec == nil {
				return ""
			}
			return rec.GetKey()
		},
		step:     func(rec *datav1.MediaRecord) int64 { return rec.GetStep() },
		precheck: (*planner.MetricPlanner).PrecheckMedia,
		build:    builder.BuildMediaRecord,
		commit: func(rec *datav1.MediaRecord) error {
			_, err := s.metrics.CommitMedia(rec, nil)
			return err
		},
	})
}

func (s *Session) upsertHistograms(ctx context.Context, records []*datav1.HistogramRecord) error {
	return commitChunks(ctx, s, records, s.prepareHistograms)
}

func (s *Session) prepareHistograms(records []*datav1.HistogramRecord) (preparedBatch, error) {
	return prepareMetricBatch(s, records, metricRecordOps[*datav1.HistogramRecord]{
		kindLabel: "histogram",
		key: func(rec *datav1.HistogramRecord) string {
			if rec == nil {
				return ""
			}
			return rec.GetKey()
		},
		step:     func(rec *datav1.HistogramRecord) int64 { return rec.GetStep() },
		precheck: (*planner.MetricPlanner).PrecheckHistogram,
		build:    builder.BuildHistogramRecord,
		commit: func(rec *datav1.HistogramRecord) error {
			_, err := s.metrics.CommitHistogram(rec, nil)
			return err
		},
	})
}

// upsertLogs 单批提交日志，超限时整批拒绝。
func (s *Session) upsertLogs(ctx context.Context, records []*terminalv1.LogRecord) error {
	if err := s.checkBatchLimit(len(records)); err != nil {
		return err
	}
	return s.commitPrepared(ctx, func() (preparedBatch, error) {
		return s.prepareLogs(records)
	})
}

func (s *Session) prepareLogs(records []*terminalv1.LogRecord) (preparedBatch, error) {
	planned := make([]*recordv1.Record, 0, len(records))
	for _, rec := range records {
		if rec == nil {
			return preparedBatch{}, fmt.Errorf("%w: log record must not be nil", ErrInvalidRecord)
		}
		planned = append(planned, builder.BuildLogRecord(s.counter, s.epoch, rec))
	}
	return preparedBatch{records: planned}, nil
}

// upsertSaves 单批提交保存记录，超限时整批拒绝。
func (s *Session) upsertSaves(ctx context.Context, records []*savev1.SaveRecord) error {
	if err := s.checkBatchLimit(len(records)); err != nil {
		return err
	}
	return s.commitPrepared(ctx, func() (preparedBatch, error) {
		return s.prepareSaves(records)
	})
}

// prepareSaves 构建保存记录：CUSTOM 类型先在本地建立文件链接，待容量预留成功后再注册 LIVE 监听。
func (s *Session) prepareSaves(records []*savev1.SaveRecord) (preparedBatch, error) {
	internal := make([]*savev1.SaveRecord, 0, len(records))
	custom := make([]*savev1.SaveRecord, 0, len(records))
	for _, rec := range records {
		if rec == nil {
			return preparedBatch{}, fmt.Errorf("%w: save record must not be nil", ErrInvalidRecord)
		}
		if rec.GetType() == savev1.SaveType_SAVE_TYPE_CUSTOM {
			custom = append(custom, rec)
			continue
		}
		internal = append(internal, rec)
	}

	planned := make([]*recordv1.Record, 0, len(records))
	for _, rec := range internal {
		built, err := builder.BuildSaveRecord(s.counter, rec)
		if err != nil {
			return preparedBatch{}, fmt.Errorf("%w: %v", ErrInvalidRecord, err)
		}
		planned = append(planned, built)
	}

	if len(custom) > 0 {
		if !s.cfg.SkipStore {
			created, err := watcher.CreateSaveLinks(custom, s.run.FilesDir())
			if err != nil {
				return preparedBatch{}, fmt.Errorf("upsert saves: create save links: %w", err)
			}
			if created > 0 {
				console.Infof("Symlinked %d files into the SwanLab run directory; "+
					"call swanlab.save again to sync new files.", created)
			}
		}
		for _, rec := range custom {
			built, err := builder.BuildSaveRecord(s.counter, rec)
			if err != nil {
				return preparedBatch{}, fmt.Errorf("%w: %v", ErrInvalidRecord, err)
			}
			planned = append(planned, built)
		}
	}
	apply := func() error {
		if s.cfg.SkipStore || len(custom) == 0 {
			return nil
		}
		// 预留成功后才注册监听，避免被拒批次留下上传副作用。
		return s.watcher.RegisterLive(custom)
	}
	return preparedBatch{records: planned, apply: apply}, nil
}
