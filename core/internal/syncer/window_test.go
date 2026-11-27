package syncer

import (
	"errors"
	"testing"

	"google.golang.org/protobuf/types/known/timestamppb"

	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
)

func ts(sec int64, nanos int32) *timestamppb.Timestamp {
	return &timestamppb.Timestamp{Seconds: sec, Nanos: nanos}
}

func i64(v int64) *int64 { return &v }

func TestFromProtoEmptyAndZeroTimestamps(t *testing.T) {
	f, err := FromProto(nil)
	if err != nil {
		t.Fatalf("FromProto(nil): %v", err)
	}
	if f.Enabled() || f.HasTime() || !f.MatchMetric(1, ts(1, 0)) {
		t.Fatalf("nil range should be disabled and match all: %+v", f)
	}

	f, err = FromProto(&corev1.SyncRange{StartTime: ts(0, 0), EndTime: ts(0, 0)})
	if err != nil {
		t.Fatalf("FromProto(zero ts): %v", err)
	}
	if f.Enabled() || f.HasTime() {
		t.Fatalf("zero timestamps should not enable filter: %+v", f)
	}
}

func TestFromProtoRejectsInvalidBounds(t *testing.T) {
	cases := []*corev1.SyncRange{
		{StartStep: i64(-1)},
		{EndStep: i64(-2)},
		{StartStep: i64(10), EndStep: i64(5)},
		{StartTime: ts(200, 0), EndTime: ts(100, 0)},
		{StartTime: ts(100, -1)},
		{EndTime: ts(100, 1_000_000_000)},
	}
	for i, tc := range cases {
		if _, err := FromProto(tc); !errors.Is(err, ErrInvalidSyncRange) {
			t.Fatalf("case %d: got %v, want ErrInvalidSyncRange", i, err)
		}
	}
}

func TestRangeFilterStepBounds(t *testing.T) {
	f, err := FromProto(&corev1.SyncRange{StartStep: i64(10), EndStep: i64(20)})
	if err != nil {
		t.Fatal(err)
	}
	if !f.Enabled() || f.HasTime() {
		t.Fatalf("unexpected flags: enabled=%v hasTime=%v", f.Enabled(), f.HasTime())
	}
	for _, tc := range []struct {
		step int64
		want bool
	}{
		{9, false},
		{10, true},
		{15, true},
		{20, true},
		{21, false},
	} {
		if got := f.MatchMetric(tc.step, nil); got != tc.want {
			t.Fatalf("MatchMetric(%d, nil) = %v, want %v", tc.step, got, tc.want)
		}
	}
	// 无时间范围时，时间维不拦截仅按时间匹配的调用方（终端日志由 HasTime() 控制）。
	if !f.MatchTime(ts(1, 0)) {
		t.Fatal("MatchTime should pass when no time bounds are set")
	}

	open, err := FromProto(&corev1.SyncRange{StartStep: i64(100)})
	if err != nil {
		t.Fatal(err)
	}
	if open.MatchMetric(99, nil) || !open.MatchMetric(100, nil) || !open.MatchMetric(10_000, nil) {
		t.Fatal("open-ended start_step mismatch")
	}
}

func TestRangeFilterTimeBoundsAndSubsecond(t *testing.T) {
	f, err := FromProto(&corev1.SyncRange{
		StartTime: ts(100, 500),
		EndTime:   ts(200, 0),
	})
	if err != nil {
		t.Fatal(err)
	}
	if !f.Enabled() || !f.HasTime() {
		t.Fatalf("expected enabled and hasTime: %+v", f)
	}
	for _, tc := range []struct {
		at   *timestamppb.Timestamp
		want bool
	}{
		{ts(100, 499), false},
		{ts(100, 500), true},
		{ts(150, 0), true},
		{ts(200, 0), true},
		{ts(200, 1), false},
	} {
		if got := f.MatchMetric(1, tc.at); got != tc.want {
			t.Fatalf("MatchMetric(1, %v) = %v, want %v", tc.at, got, tc.want)
		}
		if got := f.MatchTime(tc.at); got != tc.want {
			t.Fatalf("MatchTime(%v) = %v, want %v", tc.at, got, tc.want)
		}
	}

	// 同一微秒端点构成闭区间必须命中。
	exactTS := ts(1722470400, 1000)
	exact, err := FromProto(&corev1.SyncRange{StartTime: exactTS, EndTime: exactTS})
	if err != nil {
		t.Fatal(err)
	}
	if !exact.MatchTime(exactTS) {
		t.Fatal("exact timestamp point should match inclusive bounds")
	}
}

func TestRangeFilterCombinesStepAndTime(t *testing.T) {
	f, err := FromProto(&corev1.SyncRange{
		StartStep: i64(2),
		EndStep:   i64(3),
		StartTime: ts(100, 0),
	})
	if err != nil {
		t.Fatal(err)
	}
	if !f.MatchMetric(2, ts(100, 0)) {
		t.Fatal("step=2 ts=100 should match")
	}
	if f.MatchMetric(1, ts(100, 0)) {
		t.Fatal("step=1 should be rejected by start_step")
	}
	if f.MatchMetric(2, ts(99, 0)) {
		t.Fatal("ts=99 should be rejected by start_time")
	}
}

func TestPickTimestampFallback(t *testing.T) {
	p := ts(200, 5)
	r := ts(100, 1)
	if got := PickTimestamp(p, r); got != p {
		t.Fatalf("PickTimestamp(p, r) = %v, want payload %v", got, p)
	}
	if got := PickTimestamp(ts(0, 0), r); got != r {
		t.Fatalf("PickTimestamp(zero, r) = %v, want record %v", got, r)
	}
	if got := PickTimestamp(nil, r); got != r {
		t.Fatalf("PickTimestamp(nil, r) = %v, want record %v", got, r)
	}
	zero := PickTimestamp(nil, nil)
	if zero == nil || zero.Seconds != 0 || zero.Nanos != 0 {
		t.Fatalf("PickTimestamp(nil, nil) = %v, want zero Timestamp", zero)
	}
}
