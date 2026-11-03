package manager

import (
	"fmt"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

// summaryFromAPI 把云端摘要转换为指标恢复层使用的纯 Go 结构。
// 两侧字段一一对应，nil 输入保持 nil，交给 metrics.NewFromSummary 处理。
func summaryFromAPI(s *api.Summary) *metrics.Summary {
	if s == nil {
		return nil
	}
	return &metrics.Summary{
		Log:    summaryItemsFromAPI(s.Log),
		Media:  summaryItemsFromAPI(s.Media),
		Scalar: summaryItemsFromAPI(s.Scalar),
	}
}

func summaryItemsFromAPI(items []api.ColumnStep) []metrics.SummaryItem {
	if len(items) == 0 {
		return nil
	}
	out := make([]metrics.SummaryItem, len(items))
	for i, item := range items {
		out[i] = metrics.SummaryItem{Key: item.Key, Step: item.Step}
	}
	return out
}

// resumeModeFromProto 把 proto 枚举转换为云端接口使用的恢复策略。
// 未知取值属于协议错误，不能静默降级为 allow/never。
func resumeModeFromProto(mode runv1.ResumeMode) (api.ResumeMode, error) {
	switch mode {
	case runv1.ResumeMode_RESUME_MODE_NEVER:
		return api.ResumeNever, nil
	case runv1.ResumeMode_RESUME_MODE_ALLOW:
		return api.ResumeAllow, nil
	case runv1.ResumeMode_RESUME_MODE_MUST:
		return api.ResumeMust, nil
	default:
		return "", fmt.Errorf("unknown resume mode %d", mode)
	}
}

// runStateFromProto 把 proto 枚举转换为云端上报的终态。
// 与 Python 的 stop_experiment 一致：只有 crashed/aborted 是特殊值，
// 其余（含 RUNNING 与未知值）一律按正常结束处理。
func runStateFromProto(state runv1.RunState) api.RunState {
	switch state {
	case runv1.RunState_RUN_STATE_CRASHED:
		return api.RunStateCrashed
	case runv1.RunState_RUN_STATE_ABORTED:
		return api.RunStateAborted
	default:
		return api.RunStateFinished
	}
}
