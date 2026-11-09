// Package health 汇总站点控制面 / 数据面的探测结果并导出降级状态。
//
// 探针只判断可达性，HTTP 语义由 ProbeFunc 决定；依赖级故障（PG / Redis / Vector）
// 不在探针范围内，由写入侧失败驱动降级兜底。
package health

import "time"

// State 是聚合后的站点健康状态。
type State int32

const (
	// StateNormal 表示两个平面都未达到降级阈值，且写入侧没有未解除的降级信号。
	StateNormal State = iota
	// StateDegraded 表示至少一个平面连续不可达，或写入侧连续失败。
	StateDegraded
	// StateBackfill 预留，当前不会产生。
	StateBackfill
)

func (s State) String() string {
	switch s {
	case StateNormal:
		return "NORMAL"
	case StateDegraded:
		return "DEGRADED"
	case StateBackfill:
		return "BACKFILL"
	default:
		return "UNKNOWN"
	}
}

// Reason 是降级来源，取值与 stats 的 degraded_reason 一致。
type Reason string

const (
	// ReasonProbeControl 表示控制面（Server）连续不可达。
	ReasonProbeControl Reason = "probe:control"
	// ReasonProbeData 表示数据面（House）连续不可达。
	ReasonProbeData Reason = "probe:data"
	// ReasonWrite 表示写入侧连续可重试失败，由上层通过 Report 注入。
	ReasonWrite Reason = "write"
)

// Snapshot 是站点健康的一次快照。
type Snapshot struct {
	State  State
	Reason Reason
	// ControlReachable / DataReachable 是各平面最近一次探测结果：单次抖动会反映在这里，
	// 但不会立即降级。
	ControlReachable bool
	DataReachable    bool
}

// 默认探针节奏与阈值；SiteOptions 中的 0 值会被替换为这里的默认值。
const (
	DefaultInterval         = 5 * time.Second
	DefaultTimeout          = 2 * time.Second
	DefaultDegradeThreshold = 3
	DefaultRecoverThreshold = 2
)

// positiveOrDefault 返回 v；v <= 0 时返回默认值。
func positiveOrDefault[T ~int | ~int64](v, def T) T {
	if v <= 0 {
		return def
	}
	return v
}
