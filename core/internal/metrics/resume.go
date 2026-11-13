package metrics

import (
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
)

// Summary 是云端实验摘要中与指标恢复相关的部分。
// 纯 Go 结构：云端 JSON 的解析由 api 层完成。
type Summary struct {
	// Log 控制台日志摘要，取第一条的 Step 作为 console_epoch
	//（DeliverRunStartResponse 需要还原终端行号）。
	Log []SummaryItem
	// Media 媒体列摘要（key + 已到达的最大 step）。
	Media []SummaryItem
	// Scalar 标量列摘要；系统列（__swanlab__. 前缀）计入 global_system_step。
	Scalar []SummaryItem
	// Histogram 直方图列摘要；直方图恒为用户自定义列，计入 global_step。
	Histogram []SummaryItem
}

// SummaryItem 是单列摘要：该列在云端已记录到的最大 step。
type SummaryItem struct {
	Key  string
	Step int64
}

// ResumeResult 是 NewFromSummary 的恢复结果。
type ResumeResult struct {
	Metrics *RunMetrics
	// ConsoleEpoch 恢复的终端行号（日志摘要为空时为 0）。
	ConsoleEpoch int64
	// GlobalStep 用户侧指标（自定义列 + 媒体列）的最大 step，未恢复时为 -1。
	GlobalStep int64
	// GlobalSystemStep 系统采集侧指标的最大 step，未恢复时为 -1。
	GlobalSystemStep int64
}

// NewFromSummary 依据云端实验摘要恢复指标状态机：
//
//   - 媒体列构建 UNSPECIFIED 占位列，min_step 置为云端已到达的 step；
//   - 标量列按系统/自定义分类构建降级列，分别推进两条 step 水位线；
//   - 直方图列构建 HISTOGRAM 列，min_step 置为云端已到达的 step；
//   - 各列位图从空开始，靠 min_step 水位线拦截历史 step。
//
// 摘要中出现重复 key 视为云端数据异常，返回错误。
func NewFromSummary(s *Summary) (*ResumeResult, error) {
	return MergeSummary(nil, s)
}

// MergeSummary 把本地 Rebuild 结果与可选云端摘要合并：
//
//   - 本地已定义的列保留原列定义、位图与标量极值；若云端 minStep 更高，仅抬升其水位线；
//   - 仅存在于云端的列补充占位定义并设置 minStep；
//   - ConsoleEpoch、GlobalStep、GlobalSystemStep 取两侧最大值。
func MergeSummary(local *RebuildResult, s *Summary) (*ResumeResult, error) {
	res := &ResumeResult{Metrics: New(), GlobalStep: -1, GlobalSystemStep: -1}
	if local != nil {
		if local.Metrics != nil {
			res.Metrics = local.Metrics
		}
		res.ConsoleEpoch = local.Epoch
		res.GlobalStep = local.GlobalStep
		res.GlobalSystemStep = local.GlobalSystemStep
	}
	if s == nil {
		return res, nil
	}

	if len(s.Log) > 0 && s.Log[0].Step > res.ConsoleEpoch {
		res.ConsoleEpoch = s.Log[0].Step
	}

	seen := make(map[string]struct{}, len(s.Media)+len(s.Scalar)+len(s.Histogram))
	for _, item := range s.Media {
		if item.Step > res.GlobalStep {
			res.GlobalStep = item.Step
		}
		if err := mergeSummaryItem(res.Metrics, seen, item.Key, item.Step, KindMedia, "media", func() *columnv1.ColumnRecord {
			return builder.BuildResumeColumn(item.Key, true, false)
		}); err != nil {
			return nil, err
		}
	}

	for _, item := range s.Scalar {
		system := builder.IsSystemKey(item.Key)
		if system {
			if item.Step > res.GlobalSystemStep {
				res.GlobalSystemStep = item.Step
			}
		} else if item.Step > res.GlobalStep {
			res.GlobalStep = item.Step
		}
		if err := mergeSummaryItem(res.Metrics, seen, item.Key, item.Step, KindScalar, "scalar", func() *columnv1.ColumnRecord {
			return builder.BuildResumeColumn(item.Key, false, system)
		}); err != nil {
			return nil, err
		}
	}

	for _, item := range s.Histogram {
		if item.Step > res.GlobalStep {
			res.GlobalStep = item.Step
		}
		if err := mergeSummaryItem(res.Metrics, seen, item.Key, item.Step, KindHistogram, "histogram", func() *columnv1.ColumnRecord {
			return builder.BuildHistogramResumeColumn(item.Key)
		}); err != nil {
			return nil, err
		}
	}

	return res, nil
}

func mergeSummaryItem(
	m *RunMetrics,
	seen map[string]struct{},
	key string,
	step int64,
	kind Kind,
	label string,
	buildCol func() *columnv1.ColumnRecord,
) error {
	if _, dup := seen[key]; dup {
		return fmt.Errorf("restore %s metric %q: %w", label, key, ErrColumnExists)
	}
	seen[key] = struct{}{}
	if mt, exists := m.metrics[key]; exists {
		if step > mt.minStep {
			mt.minStep = step
		}
		return nil
	}
	if err := m.Define(key, buildCol(), step, kind); err != nil {
		return fmt.Errorf("restore %s metric %q: %w", label, key, err)
	}
	return nil
}
