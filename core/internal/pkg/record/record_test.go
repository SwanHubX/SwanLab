package record

import (
	"testing"

	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
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
