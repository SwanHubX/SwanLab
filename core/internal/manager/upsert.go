package manager

import (
	"context"
	"errors"
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

// preparedBatch 是一次两阶段提交的产物：待写入的记录与需要在预留成功后
// 才执行的指标状态变更。build 阶段只读，apply 阶段才修改 metrics。
type preparedBatch struct {
	records []*recordv1.Record
	apply   func() error
}

// UpsertColumns 定义一批列；已定义的 key 跳过，空列表按 no-op 成功返回。
func (m *Manager) UpsertColumns(ctx context.Context, sessionID string, columns []*columnv1.ColumnRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertColumns(ctx, columns)
}

// UpsertScalars 提交一批标量；新 key 自动定义列，重复/过期 step 逐条跳过。
func (m *Manager) UpsertScalars(ctx context.Context, sessionID string, records []*datav1.ScalarRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertScalars(ctx, records)
}

// UpsertMedia 提交一批媒体指标；语义与 UpsertScalars 一致。
func (m *Manager) UpsertMedia(ctx context.Context, sessionID string, records []*datav1.MediaRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertMedia(ctx, records)
}

// UpsertLogs 提交一批终端日志，epoch 从 Resume 恢复的行号继续递增。
func (m *Manager) UpsertLogs(ctx context.Context, sessionID string, records []*terminalv1.LogRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertLogs(ctx, records)
}

// UpsertSaves 提交一批保存记录：内部类型只做元数据上传，CUSTOM 额外建立
// 本地镜像并按 policy 立即上传或延迟到 Close。
func (m *Manager) UpsertSaves(ctx context.Context, sessionID string, records []*savev1.SaveRecord) error {
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.upsertSaves(ctx, records)
}

// commitPrepared 在会话锁内完成「状态校验 → 构建 → 预留 → 指标提交 → 入队」。
// build 只读，Reserve 失败时不会留下任何指标状态变更；但 build 已经分配的
// counter 编号会被消耗，与 Python 异常丢批时的跳号行为一致。
func (s *Session) commitPrepared(ctx context.Context, build func() (preparedBatch, error)) error {
	s.mu.Lock()
	defer s.mu.Unlock()

	if RunState(s.state.Load()) != RunStateRunning {
		return ErrNotRunning
	}
	batch, err := build()
	if err != nil {
		return err
	}
	if len(batch.records) == 0 {
		return nil
	}
	ticket, err := s.intake.Reserve(ctx, demandFor(batch.records, s.cfg.SkipStore))
	if err != nil {
		return upsertError(err)
	}
	if batch.apply != nil {
		if err = batch.apply(); err != nil {
			s.intake.Release(ticket)
			return fmt.Errorf("%w: %v", ErrInvalidRecord, err)
		}
	}
	if err = s.intake.Commit(ticket, batch.records); err != nil {
		return upsertError(err)
	}
	return nil
}

// upsertError 把 intake 的容量错误映射为 manager 的领域错误。
func upsertError(err error) error {
	switch {
	case errors.Is(err, buffer.ErrOverloaded), errors.Is(err, buffer.ErrBatchTooLarge):
		return fmt.Errorf("%w: %v", ErrOverloaded, err)
	case errors.Is(err, buffer.ErrClosed):
		return fmt.Errorf("%w: %w", ErrNotRunning, err)
	default:
		return err
	}
}

func (s *Session) upsertColumns(ctx context.Context, columns []*columnv1.ColumnRecord) error {
	size := s.chunkSize()
	for start := 0; start < len(columns); start += size {
		chunk := columns[start:min(start+size, len(columns))]
		if err := s.commitPrepared(ctx, func() (preparedBatch, error) {
			return s.prepareColumns(chunk)
		}); err != nil {
			return err
		}
	}
	return nil
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
		defs = append(defs, columnPlan{key: key, column: column, scalar: isScalarColumn(column)})
		records = append(records, builder.BuildColumnRecord(s.counter, column))
	}
	apply := func() error {
		for _, def := range defs {
			var err error
			if def.scalar {
				err = s.metrics.DefineScalar(def.key, def.column, -1)
			} else {
				err = s.metrics.DefineMedia(def.key, def.column, -1)
			}
			if err != nil {
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

func isScalarColumn(column *columnv1.ColumnRecord) bool {
	return column.GetColumnType() == columnv1.ColumnType_COLUMN_TYPE_SCALAR
}

func (s *Session) upsertScalars(ctx context.Context, records []*datav1.ScalarRecord) error {
	size := s.chunkSize()
	for start := 0; start < len(records); start += size {
		chunk := records[start:min(start+size, len(records))]
		if err := s.commitPrepared(ctx, func() (preparedBatch, error) {
			return s.prepareScalars(chunk)
		}); err != nil {
			return err
		}
	}
	return nil
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
	size := s.chunkSize()
	for start := 0; start < len(records); start += size {
		chunk := records[start:min(start+size, len(records))]
		if err := s.commitPrepared(ctx, func() (preparedBatch, error) {
			return s.prepareMedia(chunk)
		}); err != nil {
			return err
		}
	}
	return nil
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

func (s *Session) upsertLogs(ctx context.Context, records []*terminalv1.LogRecord) error {
	size := s.chunkSize()
	for start := 0; start < len(records); start += size {
		chunk := records[start:min(start+size, len(records))]
		if err := s.commitPrepared(ctx, func() (preparedBatch, error) {
			return s.prepareLogs(chunk)
		}); err != nil {
			return err
		}
	}
	return nil
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

func (s *Session) upsertSaves(ctx context.Context, records []*savev1.SaveRecord) error {
	size := s.chunkSize()
	for start := 0; start < len(records); start += size {
		chunk := records[start:min(start+size, len(records))]
		if err := s.commitPrepared(ctx, func() (preparedBatch, error) {
			return s.prepareSaves(chunk)
		}); err != nil {
			return err
		}
	}
	return nil
}

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
			if err = s.watcher.RegisterLive(custom); err != nil {
				return preparedBatch{}, fmt.Errorf("upsert saves: watch live saves: %w", err)
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
	return preparedBatch{records: planned}, nil
}

// demandFor 计算一批记录对本地与上传两条通道的需求。
// 只有 CUSTOM + SAVE_POLICY_END 会被延迟到 Close，与 buffer.routeRecords 同源。
func demandFor(records []*recordv1.Record, skipStore bool) buffer.Demand {
	var demand buffer.Demand
	if !skipStore {
		demand.Store = len(records)
	}
	for _, rec := range records {
		size := record.Size(rec)
		switch record.TypeOf(rec) {
		case record.TypeStart, record.TypeFinish:
		case record.TypeSave:
			if deferredSave(rec) {
				demand.Deferred++
				demand.DeferredBytes += size
				continue
			}
			demand.Upload++
			demand.Bytes += size
		default:
			demand.Upload++
			demand.Bytes += size
		}
	}
	return demand
}

func deferredSave(rec *recordv1.Record) bool {
	save := rec.GetSave()
	return save.GetPolicy() == savev1.SavePolicy_SAVE_POLICY_END &&
		save.GetType() == savev1.SaveType_SAVE_TYPE_CUSTOM
}

func (s *Session) chunkSize() int {
	if s.cfg.RecordBatch <= 0 {
		return 1
	}
	return int(s.cfg.RecordBatch)
}

// columnPlan 记录本批次内新定义的列；apply 阶段按该顺序落地，保证与记录顺序一致。
type columnPlan struct {
	key    string
	column *columnv1.ColumnRecord
	scalar bool
}

// metricPlanner 在 Reserve 之前对一批标量/媒体做只读预检，并在本地模拟同批内
// 已定义列与新 step，避免同一 key 的多条记录重复产出列或重复接收 step。
type metricPlanner struct {
	metrics     *metrics.RunMetrics
	sectionRule int32
	columns     map[string]*columnv1.ColumnRecord
	plans       []columnPlan
	steps       map[string]map[int64]struct{}
}

func newMetricPlanner(m *metrics.RunMetrics, sectionRule int32) *metricPlanner {
	return &metricPlanner{
		metrics:     m,
		sectionRule: sectionRule,
		columns:     make(map[string]*columnv1.ColumnRecord),
		steps:       make(map[string]map[int64]struct{}),
	}
}

func (p *metricPlanner) precheckScalar(rec *datav1.ScalarRecord) metrics.Precheck {
	if column, ok := p.columns[rec.GetKey()]; ok {
		return p.remember(rec.GetKey(), column.GetColumnType(), rec.GetStep(), true,
			pendingVerdict(rec.GetKey(), column.GetColumnType(), rec.GetType(), rec.GetStep()))
	}
	return p.remember(rec.GetKey(), rec.GetType(), rec.GetStep(), true, p.metrics.PrecheckScalar(rec))
}

func (p *metricPlanner) precheckMedia(rec *datav1.MediaRecord) metrics.Precheck {
	if column, ok := p.columns[rec.GetKey()]; ok {
		return p.remember(rec.GetKey(), column.GetColumnType(), rec.GetStep(), false,
			pendingVerdict(rec.GetKey(), column.GetColumnType(), rec.GetType(), rec.GetStep()))
	}
	return p.remember(rec.GetKey(), rec.GetType(), rec.GetStep(), false, p.metrics.PrecheckMedia(rec))
}

// remember 登记本批次内新建的列与已接纳的 step，并把同批重复 step 改判为跳过。
func (p *metricPlanner) remember(
	key string,
	typ columnv1.ColumnType,
	step int64,
	scalar bool,
	verdict metrics.Precheck,
) metrics.Precheck {
	if verdict.NeedsColumn {
		column := builder.BuildAutoColumn(p.sectionRule, key, typ)
		p.columns[key] = column
		// 定义路径由 Upsert RPC 决定，与 Python define_scalar/define_media 对齐；
		// 不能从 column 类型反推（异常输入下二者可能不一致）。
		p.plans = append(p.plans, columnPlan{key: key, column: column, scalar: scalar})
		p.steps[key] = make(map[int64]struct{})
	}
	if verdict.Err != nil || verdict.Decision != metrics.DecisionAccept {
		return verdict
	}
	if p.seenStep(key, step) {
		verdict.Decision = metrics.DecisionSkip
		return verdict
	}
	p.markStep(key, step)
	return verdict
}

func (p *metricPlanner) seenStep(key string, step int64) bool {
	_, ok := p.steps[key][step]
	return ok
}

func (p *metricPlanner) markStep(key string, step int64) {
	if p.steps[key] == nil {
		p.steps[key] = make(map[int64]struct{})
	}
	p.steps[key][step] = struct{}{}
}

// defineAll 把本批次内自动定义的列写入 metrics；必须在 Reserve 成功之后调用。
func (p *metricPlanner) defineAll() error {
	for _, plan := range p.plans {
		var err error
		if plan.scalar {
			err = p.metrics.DefineScalar(plan.key, plan.column, -1)
		} else {
			err = p.metrics.DefineMedia(plan.key, plan.column, -1)
		}
		if err != nil {
			return err
		}
	}
	return nil
}

// pendingVerdict 对「本批次内刚定义、尚未写入 metrics」的 key 复刻
// freshVerdict 与类型校验，保证同批重复 step / 类型冲突的结论与逐条提交一致。
func pendingVerdict(
	key string,
	columnType columnv1.ColumnType,
	recordType columnv1.ColumnType,
	step int64,
) metrics.Precheck {
	if columnType != recordType {
		return metrics.Precheck{
			Decision: metrics.DecisionSkip,
			Err: fmt.Errorf("%w: metric '%s' has already been defined as %s, not %s",
				metrics.ErrTypeMismatch, key, columnType, recordType),
		}
	}
	switch {
	case step <= -1:
		return metrics.Precheck{Decision: metrics.DecisionSkip}
	case step > metrics.MaxStep:
		return metrics.Precheck{
			Decision: metrics.DecisionSkip,
			Err:      fmt.Errorf("step %d: %w (max %d)", step, metrics.ErrStepLimit, metrics.MaxStep),
		}
	case step < 0:
		return metrics.Precheck{
			Decision: metrics.DecisionSkip,
			Err:      fmt.Errorf("step %d: %w", step, metrics.ErrStepNegative),
		}
	}
	return metrics.Precheck{Decision: metrics.DecisionAccept}
}
