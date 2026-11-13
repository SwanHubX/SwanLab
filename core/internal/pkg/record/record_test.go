package record

import (
	"testing"

	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func TestTypeOf(t *testing.T) {
	tests := []struct {
		name string
		in   *recordv1.Record
		want Type
	}{
		{"nil", nil, ""},
		{"unset", &recordv1.Record{}, ""},
		{"start", &recordv1.Record{RecordType: &recordv1.Record_Start{}}, TypeStart},
		{"finish", &recordv1.Record{RecordType: &recordv1.Record_Finish{}}, TypeFinish},
		{"column", &recordv1.Record{RecordType: &recordv1.Record_Column{}}, TypeColumn},
		{"scalar", &recordv1.Record{RecordType: &recordv1.Record_Scalar{}}, TypeScalar},
		{"media", &recordv1.Record{RecordType: &recordv1.Record_Media{}}, TypeMedia},
		{"histogram", &recordv1.Record{RecordType: &recordv1.Record_Histogram{}}, TypeHistogram},
		{"log", &recordv1.Record{RecordType: &recordv1.Record_Log{}}, TypeLog},
		{"save", &recordv1.Record{RecordType: &recordv1.Record_Save{}}, TypeSave},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := TypeOf(tt.in); got != tt.want {
				t.Fatalf("TypeOf() = %q, want %q", got, tt.want)
			}
		})
	}
}

func TestSize(t *testing.T) {
	if got := Size(nil); got != 0 {
		t.Fatalf("Size(nil) = %d, want 0", got)
	}
	if got := Size(&recordv1.Record{Num: 1}); got <= 0 {
		t.Fatalf("Size() = %d, want > 0", got)
	}
}

func TestUploadRoutingPredicates(t *testing.T) {
	customEnd := &recordv1.Record{
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{Policy: savev1.SavePolicy_SAVE_POLICY_END, Type: savev1.SaveType_SAVE_TYPE_CUSTOM},
		},
	}
	configEnd := &recordv1.Record{
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{Policy: savev1.SavePolicy_SAVE_POLICY_END, Type: savev1.SaveType_SAVE_TYPE_CONFIG},
		},
	}
	customLive := &recordv1.Record{
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{Policy: savev1.SavePolicy_SAVE_POLICY_LIVE, Type: savev1.SaveType_SAVE_TYPE_CUSTOM},
		},
	}

	cases := []struct {
		name      string
		rec       *recordv1.Record
		deferred  bool
		immediate bool
	}{
		{"nil", nil, false, false},
		{"start", &recordv1.Record{RecordType: &recordv1.Record_Start{}}, false, false},
		{"finish", &recordv1.Record{RecordType: &recordv1.Record_Finish{}}, false, false},
		{"scalar", &recordv1.Record{RecordType: &recordv1.Record_Scalar{}}, false, true},
		{"custom_end", customEnd, true, false},
		{"config_end", configEnd, false, true},
		{"custom_live", customLive, false, true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := IsDeferredSave(tc.rec); got != tc.deferred {
				t.Fatalf("IsDeferredSave = %v, want %v", got, tc.deferred)
			}
			if got := IsImmediateUploadable(tc.rec); got != tc.immediate {
				t.Fatalf("IsImmediateUploadable = %v, want %v", got, tc.immediate)
			}
		})
	}
}
