package metrics

import (
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
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
//   - 各列位图从空开始，靠 min_step 水位线拦截历史 step。
//
// 摘要中出现重复 key 视为云端数据异常，返回错误。
func NewFromSummary(s *Summary) (*ResumeResult, error) {
	res := &ResumeResult{Metrics: New(), GlobalStep: -1, GlobalSystemStep: -1}
	if s == nil {
		return res, nil
	}

	if len(s.Log) > 0 {
		res.ConsoleEpoch = s.Log[0].Step
	}

	for _, item := range s.Media {
		if item.Step > res.GlobalStep {
			res.GlobalStep = item.Step
		}
		col := builder.BuildResumeColumn(item.Key, true, false)
		if err := res.Metrics.DefineMedia(item.Key, col, item.Step); err != nil {
			return nil, fmt.Errorf("restore media metric %q: %w", item.Key, err)
		}
	}

	for _, item := range s.Scalar {
		system := builder.IsSystemKey(item.Key)
		if system {
			if item.Step > res.GlobalSystemStep {
				res.GlobalSystemStep = item.Step
			}
		} else {
			if item.Step > res.GlobalStep {
				res.GlobalStep = item.Step
			}
		}
		col := builder.BuildResumeColumn(item.Key, false, system)
		if err := res.Metrics.DefineScalar(item.Key, col, item.Step); err != nil {
			return nil, fmt.Errorf("restore scalar metric %q: %w", item.Key, err)
		}
	}

	return res, nil
}
