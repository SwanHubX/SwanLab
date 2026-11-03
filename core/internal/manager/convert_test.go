package manager

import (
	"reflect"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func TestSummaryFromAPI(t *testing.T) {
	got := summaryFromAPI(&api.Summary{
		Log:    []api.ColumnStep{{Key: "console", Step: 3}},
		Media:  []api.ColumnStep{{Key: "image", Step: 7}},
		Scalar: []api.ColumnStep{{Key: "loss", Step: 9}, {Key: "__swanlab__.gpu", Step: 11}},
	})
	want := &metrics.Summary{
		Log:    []metrics.SummaryItem{{Key: "console", Step: 3}},
		Media:  []metrics.SummaryItem{{Key: "image", Step: 7}},
		Scalar: []metrics.SummaryItem{{Key: "loss", Step: 9}, {Key: "__swanlab__.gpu", Step: 11}},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("summary = %#v, want %#v", got, want)
	}
}

func TestSummaryFromAPINilAndEmpty(t *testing.T) {
	if got := summaryFromAPI(nil); got != nil {
		t.Fatalf("nil summary = %#v, want nil", got)
	}
	if got := summaryFromAPI(&api.Summary{}); !reflect.DeepEqual(got, &metrics.Summary{}) {
		t.Fatalf("empty summary = %#v, want empty", got)
	}
}

func TestResumeModeFromProto(t *testing.T) {
	cases := []struct {
		in   runv1.ResumeMode
		want api.ResumeMode
	}{
		{runv1.ResumeMode_RESUME_MODE_NEVER, api.ResumeNever},
		{runv1.ResumeMode_RESUME_MODE_ALLOW, api.ResumeAllow},
		{runv1.ResumeMode_RESUME_MODE_MUST, api.ResumeMust},
	}
	for _, tc := range cases {
		got, err := resumeModeFromProto(tc.in)
		if err != nil {
			t.Fatalf("resumeModeFromProto(%v): %v", tc.in, err)
		}
		if got != tc.want {
			t.Fatalf("resumeModeFromProto(%v) = %q, want %q", tc.in, got, tc.want)
		}
	}
	if _, err := resumeModeFromProto(runv1.ResumeMode(99)); err == nil {
		t.Fatal("unknown resume mode must return an error")
	}
}

func TestRunStateFromProto(t *testing.T) {
	cases := []struct {
		in   runv1.RunState
		want api.RunState
	}{
		// RUNNING 不是合法的收尾状态，Python 侧同样按正常结束上报。
		{runv1.RunState_RUN_STATE_RUNNING, api.RunStateFinished},
		{runv1.RunState_RUN_STATE_FINISHED, api.RunStateFinished},
		{runv1.RunState_RUN_STATE_CRASHED, api.RunStateCrashed},
		{runv1.RunState_RUN_STATE_ABORTED, api.RunStateAborted},
		{runv1.RunState(99), api.RunStateFinished},
	}
	for _, tc := range cases {
		if got := runStateFromProto(tc.in); got != tc.want {
			t.Fatalf("runStateFromProto(%v) = %q, want %q", tc.in, got, tc.want)
		}
	}
}
