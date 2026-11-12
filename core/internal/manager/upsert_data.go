// 本文件承载 column/scalar/media/histogram/log/save 记录的预检与构建。
package manager

import (
	"context"
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
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

func (s *Session) prepareColumns(columns []*columnv1.ColumnRecord) (preparedBatch, error) {
	records := make([]*recordv1.Record, 0, len(columns))
	defs := make([]columnPlan, 0, len(columns))
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
		defs = append(defs, columnPlan{key: key, column: column, kind: columnKind(column)})
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

func (s *Session) upsertScalars(ctx context.Context, records []*datav1.ScalarRecord) error {
	return commitChunks(ctx, s, records, s.prepareScalars)
}

func (s *Session) prepareScalars(records []*datav1.ScalarRecord) (preparedBatch, error) {
	planner := newMetricPlanner(s.metrics, s.cfg.SectionRule)
	planned := make([]*recordv1.Record, 0, len(records))
	accepted := make([]*datav1.ScalarRecord, 0, len(records))
	for _, rec := range records {
		if rec == nil || rec.GetKey() == "" {
			return preparedBatch{}, fmt.Errorf("%w: scalar key is required", ErrInvalidRecord)
		}
		verdict := planner.precheckScalar(rec)
		if verdict.NeedsColumn {
			planned = append(planned, builder.BuildColumnRecord(s.counter, planner.columns[rec.GetKey()]))
		}
		switch {
		case verdict.Err != nil:
			console.Warningf("metric %q step %d skipped: %v", rec.GetKey(), rec.GetStep(), verdict.Err)
		case verdict.Decision == metrics.DecisionAccept:
			planned = append(planned, builder.BuildScalarRecord(s.counter, rec))
			accepted = append(accepted, rec)
		default:
			console.Debugf("Metric '%s' at step %d was skipped because it is duplicate or too old.",
				rec.GetKey(), rec.GetStep())
		}
	}
	apply := func() error {
		if err := planner.defineAll(); err != nil {
			return err
		}
		for _, rec := range accepted {
			if _, err := s.metrics.CommitScalar(rec, nil); err != nil {
				return err
			}
		}
		return nil
	}
	return preparedBatch{records: planned, apply: apply}, nil
}

func (s *Session) upsertMedia(ctx context.Context, records []*datav1.MediaRecord) error {
	return commitChunks(ctx, s, records, s.prepareMedia)
}

func (s *Session) prepareMedia(records []*datav1.MediaRecord) (preparedBatch, error) {
	planner := newMetricPlanner(s.metrics, s.cfg.SectionRule)
	planned := make([]*recordv1.Record, 0, len(records))
	accepted := make([]*datav1.MediaRecord, 0, len(records))
	for _, rec := range records {
		if rec == nil || rec.GetKey() == "" {
			return preparedBatch{}, fmt.Errorf("%w: media key is required", ErrInvalidRecord)
		}
		verdict := planner.precheckMedia(rec)
		if verdict.NeedsColumn {
			planned = append(planned, builder.BuildColumnRecord(s.counter, planner.columns[rec.GetKey()]))
		}
		switch {
		case verdict.Err != nil:
			console.Warningf("media metric %q step %d skipped: %v", rec.GetKey(), rec.GetStep(), verdict.Err)
		case verdict.Decision == metrics.DecisionAccept:
			planned = append(planned, builder.BuildMediaRecord(s.counter, rec))
			accepted = append(accepted, rec)
		default:
			console.Debugf("Media metric '%s' at step %d was skipped because it is duplicate or too old.",
				rec.GetKey(), rec.GetStep())
		}
	}
	apply := func() error {
		if err := planner.defineAll(); err != nil {
			return err
		}
		for _, rec := range accepted {
			if _, err := s.metrics.CommitMedia(rec, nil); err != nil {
				return err
			}
		}
		return nil
	}
	return preparedBatch{records: planned, apply: apply}, nil
}

func (s *Session) upsertHistograms(ctx context.Context, records []*datav1.HistogramRecord) error {
	return commitChunks(ctx, s, records, s.prepareHistograms)
}

func (s *Session) prepareHistograms(records []*datav1.HistogramRecord) (preparedBatch, error) {
	planner := newMetricPlanner(s.metrics, s.cfg.SectionRule)
	planned := make([]*recordv1.Record, 0, len(records))
	accepted := make([]*datav1.HistogramRecord, 0, len(records))
	for _, rec := range records {
		if rec == nil || rec.GetKey() == "" {
			return preparedBatch{}, fmt.Errorf("%w: histogram key is required", ErrInvalidRecord)
		}
		verdict := planner.precheckHistogram(rec)
		if verdict.NeedsColumn {
			planned = append(planned, builder.BuildColumnRecord(s.counter, planner.columns[rec.GetKey()]))
		}
		switch {
		case verdict.Err != nil:
			console.Warningf("histogram metric %q step %d skipped: %v", rec.GetKey(), rec.GetStep(), verdict.Err)
		case verdict.Decision == metrics.DecisionAccept:
			planned = append(planned, builder.BuildHistogramRecord(s.counter, rec))
			accepted = append(accepted, rec)
		default:
			console.Debugf("Histogram metric '%s' at step %d was skipped because it is duplicate or too old.",
				rec.GetKey(), rec.GetStep())
		}
	}
	apply := func() error {
		if err := planner.defineAll(); err != nil {
			return err
		}
		for _, rec := range accepted {
			if _, err := s.metrics.CommitHistogram(rec, nil); err != nil {
				return err
			}
		}
		return nil
	}
	return preparedBatch{records: planned, apply: apply}, nil
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
