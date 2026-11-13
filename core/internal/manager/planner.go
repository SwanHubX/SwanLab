// 本文件承载批内指标规划与通道容量需求计算。
package manager

import (
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// demandFor 计算一批记录对本地与上传两条通道的需求。
// local/offline 只写本地；online 中只有 CUSTOM + SAVE_POLICY_END 会延迟到 Close，
// 与 buffer.routeRecords 同源。
func demandFor(records []*recordv1.Record, cfg config.Settings) buffer.Demand {
	var demand buffer.Demand
	if !cfg.SkipStore {
		demand.Store = len(records)
	}
	if cfg.NormalizedMode() != config.ModeOnline {
		return demand
	}
	for _, rec := range records {
		switch {
		case record.IsDeferredSave(rec):
			demand.Deferred++
		case record.IsImmediateUploadable(rec):
			demand.Upload++
			demand.Bytes += record.Size(rec)
		}
	}
	return demand
}

// columnPlan 记录本批次内新定义的列；apply 阶段按该顺序落地，保证与记录顺序一致。
type columnPlan struct {
	key    string
	column *columnv1.ColumnRecord
	kind   metrics.Kind
}

// metricPlanner 在 Reserve 之前对一批标量/媒体/直方图做只读预检，并在本地模拟同批内
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
		return p.remember(rec.GetKey(), column.GetColumnType(), rec.GetStep(), metrics.KindScalar,
			pendingVerdict(rec.GetKey(), column.GetColumnType(), rec.GetType(), rec.GetStep()))
	}
	return p.remember(rec.GetKey(), rec.GetType(), rec.GetStep(), metrics.KindScalar, p.metrics.PrecheckScalar(rec))
}

func (p *metricPlanner) precheckMedia(rec *datav1.MediaRecord) metrics.Precheck {
	if column, ok := p.columns[rec.GetKey()]; ok {
		return p.remember(rec.GetKey(), column.GetColumnType(), rec.GetStep(), metrics.KindMedia,
			pendingVerdict(rec.GetKey(), column.GetColumnType(), rec.GetType(), rec.GetStep()))
	}
	return p.remember(rec.GetKey(), rec.GetType(), rec.GetStep(), metrics.KindMedia, p.metrics.PrecheckMedia(rec))
}

func (p *metricPlanner) precheckHistogram(rec *datav1.HistogramRecord) metrics.Precheck {
	if column, ok := p.columns[rec.GetKey()]; ok {
		return p.remember(rec.GetKey(), column.GetColumnType(), rec.GetStep(), metrics.KindHistogram,
			pendingVerdict(rec.GetKey(), column.GetColumnType(), rec.GetType(), rec.GetStep()))
	}
	return p.remember(rec.GetKey(), rec.GetType(), rec.GetStep(), metrics.KindHistogram,
		p.metrics.PrecheckHistogram(rec))
}

// remember 登记本批次内新建的列与已接纳的 step，并把同批重复 step 改判为跳过。
func (p *metricPlanner) remember(
	key string,
	typ columnv1.ColumnType,
	step int64,
	kind metrics.Kind,
	verdict metrics.Precheck,
) metrics.Precheck {
	if verdict.NeedsColumn {
		column := builder.BuildAutoColumn(p.sectionRule, key, typ)
		p.columns[key] = column
		// 定义路径由 Upsert RPC 决定，与 Python define_scalar/define_media/define_histogram 对齐；
		// 不能从 column 类型反推（异常输入下二者可能不一致）。
		p.plans = append(p.plans, columnPlan{key: key, column: column, kind: kind})
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
		if err := p.metrics.Define(plan.key, plan.column, -1, plan.kind); err != nil {
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
