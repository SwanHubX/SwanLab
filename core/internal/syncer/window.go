// Package syncer 实现离线 run-*.swanlab 文件的云端同步编排与范围过滤。
package syncer

import (
	"errors"
	"fmt"
	"time"

	"google.golang.org/protobuf/types/known/timestamppb"

	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
)

// ErrInvalidSyncRange 表示 SyncRange 参数越界、上下界倒置或时间戳非法。
var ErrInvalidSyncRange = errors.New("invalid sync range")

// RangeFilter 按 step 范围与时间范围过滤回放记录，任一维度未设置表示不限制。
//
// step 范围仅对 scalar/media/histogram 指标记录生效；时间范围对指标与终端日志生效，
// 两维度同时设置时取交集，所有上下界均为闭区间。
type RangeFilter struct {
	hasStartStep bool
	startStep    int64
	hasEndStep   bool
	endStep      int64

	hasStartTime bool
	startTime    time.Time
	hasEndTime   bool
	endTime      time.Time
}

// FromProto 从 SyncRange 协议对象构造过滤器；nil 表示全量同步。
// 时间戳为 nil 或 0（epoch）时视为未设置，避免把未赋值字段当作 1970 年过滤。
func FromProto(r *corev1.SyncRange) (RangeFilter, error) {
	if r == nil {
		return RangeFilter{}, nil
	}
	var f RangeFilter
	if r.StartStep != nil {
		if *r.StartStep < 0 {
			return RangeFilter{}, fmt.Errorf("%w: start_step must be non-negative (%d)", ErrInvalidSyncRange, *r.StartStep)
		}
		f.hasStartStep = true
		f.startStep = *r.StartStep
	}
	if r.EndStep != nil {
		if *r.EndStep < 0 {
			return RangeFilter{}, fmt.Errorf("%w: end_step must be non-negative (%d)", ErrInvalidSyncRange, *r.EndStep)
		}
		f.hasEndStep = true
		f.endStep = *r.EndStep
	}
	if f.hasStartStep && f.hasEndStep && f.startStep > f.endStep {
		return RangeFilter{}, fmt.Errorf("%w: start_step (%d) must be <= end_step (%d)", ErrInvalidSyncRange, f.startStep, f.endStep)
	}
	if isSetTimestamp(r.StartTime) {
		if err := r.StartTime.CheckValid(); err != nil {
			return RangeFilter{}, fmt.Errorf("%w: invalid start_time: %v", ErrInvalidSyncRange, err)
		}
		f.hasStartTime = true
		f.startTime = r.StartTime.AsTime()
	}
	if isSetTimestamp(r.EndTime) {
		if err := r.EndTime.CheckValid(); err != nil {
			return RangeFilter{}, fmt.Errorf("%w: invalid end_time: %v", ErrInvalidSyncRange, err)
		}
		f.hasEndTime = true
		f.endTime = r.EndTime.AsTime()
	}
	if f.hasStartTime && f.hasEndTime && f.startTime.After(f.endTime) {
		return RangeFilter{}, fmt.Errorf("%w: start_time must be <= end_time", ErrInvalidSyncRange)
	}
	return f, nil
}

// Enabled 报告是否启用了任意维度的范围过滤。
func (f RangeFilter) Enabled() bool {
	return f.hasStartStep || f.hasEndStep || f.hasStartTime || f.hasEndTime
}

// HasTime 报告是否设置了时间范围。终端日志仅在设置时间范围时参与回放。
func (f RangeFilter) HasTime() bool {
	return f.hasStartTime || f.hasEndTime
}

// MatchMetric 判断携带 step 的指标记录是否落在范围闭区间内。
func (f RangeFilter) MatchMetric(step int64, ts *timestamppb.Timestamp) bool {
	if f.hasStartStep && step < f.startStep {
		return false
	}
	if f.hasEndStep && step > f.endStep {
		return false
	}
	return f.MatchTime(ts)
}

// MatchTime 判断时间戳是否落在时间闭区间内；ts 为 nil 时不做时间拦截。
func (f RangeFilter) MatchTime(ts *timestamppb.Timestamp) bool {
	if !f.HasTime() || ts == nil {
		return true
	}
	t := ts.AsTime()
	if f.hasStartTime && t.Before(f.startTime) {
		return false
	}
	if f.hasEndTime && t.After(f.endTime) {
		return false
	}
	return true
}

// PickTimestamp 优先选用业务载荷时间，未设置（nil 或零值）时回退到记录落盘时间。
func PickTimestamp(payload, rec *timestamppb.Timestamp) *timestamppb.Timestamp {
	if isSetTimestamp(payload) {
		return payload
	}
	if rec != nil {
		return rec
	}
	return &timestamppb.Timestamp{}
}

func isSetTimestamp(ts *timestamppb.Timestamp) bool {
	return ts != nil && (ts.Seconds != 0 || ts.Nanos != 0)
}
