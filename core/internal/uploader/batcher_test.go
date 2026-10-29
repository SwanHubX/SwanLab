package uploader

import (
	"context"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func newScalar(num int64, key string) *recordv1.Record {
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Scalar{Scalar: &datav1.ScalarRecord{Key: key}},
	}
}

func newColumn(num int64, key string) *recordv1.Record {
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Column{Column: &columnv1.ColumnRecord{ColumnKey: key}},
	}
}

func newMedia(num int64, key string) *recordv1.Record {
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Media{Media: &datav1.MediaRecord{Key: key}},
	}
}

func newLog(num int64, line string) *recordv1.Record {
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Log{Log: &terminalv1.LogRecord{Line: line}},
	}
}

func newSave(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Save{Save: &savev1.SaveRecord{}},
	}
}

func newStart(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}},
	}
}

func newFinish(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num:        num,
		RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{}},
	}
}

// newTestBuffer 模拟 intake 已为每条记录预留配额，再构造缓冲区。
func newTestBuffer(t *testing.T, owned ...*recordv1.Record) (*recordBuffer, *budget.Budget, *Tracker) {
	t.Helper()
	ledger := budget.New(budget.Limits{})
	for _, rec := range owned {
		if err := ledger.Acquire(context.Background(), 1, record.Size(rec)); err != nil {
			t.Fatalf("acquire: %v", err)
		}
	}
	tracker := NewTracker()
	return newRecordBuffer(ledger, tracker), ledger, tracker
}

func TestRecordBufferAddAndDrain(t *testing.T) {
	first, second := newScalar(1, "a"), newScalar(2, "b")
	buf, _, tracker := newTestBuffer(t, first, second)

	if added, replaced := buf.add(first); !added || replaced {
		t.Fatalf("first add = (%v, %v), want (true, false)", added, replaced)
	}
	if added, replaced := buf.add(second); !added || replaced {
		t.Fatalf("second add = (%v, %v), want (true, false)", added, replaced)
	}
	if got := buf.len(); got != 2 {
		t.Fatalf("len = %d, want 2", got)
	}
	if got := tracker.Snapshot().TotalRecords; got != 2 {
		t.Fatalf("total records = %d, want 2", got)
	}

	drained := buf.drain()
	if len(drained) != 2 || drained[0] != first || drained[1] != second {
		t.Fatalf("drain = %v, want insertion order", drained)
	}
	if got := buf.len(); got != 0 {
		t.Fatalf("len after drain = %d, want 0", got)
	}
	if again := buf.drain(); again != nil {
		t.Fatalf("second drain = %v, want nil", again)
	}
}

func TestRecordBufferOverwriteKeepsPosition(t *testing.T) {
	first, second, third := newScalar(1, "a"), newScalar(2, "b"), newScalar(3, "c")
	replacement := newScalar(2, "b2")
	buf, _, tracker := newTestBuffer(t, first, second, third, replacement)

	buf.add(first)
	buf.add(second)
	buf.add(third)

	if added, replaced := buf.add(replacement); added || !replaced {
		t.Fatalf("overwrite = (%v, %v), want (false, true)", added, replaced)
	}
	if got := buf.len(); got != 3 {
		t.Fatalf("len after overwrite = %d, want 3", got)
	}

	drained := buf.drain()
	want := []*recordv1.Record{first, replacement, third}
	for i, rec := range want {
		if drained[i] != rec {
			t.Fatalf("drained[%d] = %v, want %v", i, drained[i], rec)
		}
	}
	if got := tracker.Snapshot().TotalRecords; got != 3 {
		t.Fatalf("total records = %d, want 3 (overwrite must not count)", got)
	}
}

func TestRecordBufferOverwriteReleasesBudget(t *testing.T) {
	old := newScalar(1, "old")
	replacement := newScalar(1, "a-much-longer-replacement-payload")
	other := newScalar(2, "other")
	if record.Size(old) == record.Size(replacement) {
		t.Fatal("test records must differ in size")
	}
	buf, ledger, _ := newTestBuffer(t, old, other, replacement)

	buf.add(old)
	buf.add(other)
	if used, _ := ledger.Usage(); used != 3 {
		t.Fatalf("usage before overwrite = %d records, want 3", used)
	}

	buf.add(replacement)
	used, bytes := ledger.Usage()
	if used != 2 {
		t.Fatalf("usage after overwrite = %d records, want 2", used)
	}
	if want := record.Size(other) + record.Size(replacement); bytes != want {
		t.Fatalf("usage after overwrite = %d bytes, want %d", bytes, want)
	}
}

func TestGroupRecordsOrder(t *testing.T) {
	records := []*recordv1.Record{
		newScalar(2, "s1"),
		newLog(3, "l1"),
		newColumn(1, "c1"),
		newScalar(4, "s2"),
		newSave(5),
		newMedia(6, "m1"),
		newScalar(7, "s3"),
	}

	groups := groupRecords(records)
	wantTypes := []record.Type{
		record.TypeColumn,
		record.TypeScalar,
		record.TypeMedia,
		record.TypeLog,
		record.TypeSave,
	}
	if len(groups) != len(wantTypes) {
		t.Fatalf("groups = %d, want %d", len(groups), len(wantTypes))
	}
	for i, typ := range wantTypes {
		if groups[i].typ != typ {
			t.Fatalf("group[%d] = %s, want %s", i, groups[i].typ, typ)
		}
	}
	if len(groups[0].records) != 1 || groups[0].records[0] != records[2] {
		t.Fatalf("column group = %v, want [c1]", groups[0].records)
	}

	scalars := groups[1].records
	if len(scalars) != 3 || scalars[0] != records[0] || scalars[1] != records[3] || scalars[2] != records[6] {
		t.Fatalf("scalar group order = %v, want s1,s2,s3", scalars)
	}
}

func TestGroupRecordsUnknownTypesTrailing(t *testing.T) {
	start1, start2, finish := newStart(10), newStart(12), newFinish(11)
	groups := groupRecords([]*recordv1.Record{start1, newScalar(1, "a"), finish, start2})

	wantTypes := []record.Type{record.TypeScalar, record.TypeStart, record.TypeFinish}
	if len(groups) != len(wantTypes) {
		t.Fatalf("groups = %d, want %d", len(groups), len(wantTypes))
	}
	for i, typ := range wantTypes {
		if groups[i].typ != typ {
			t.Fatalf("group[%d] = %s, want %s", i, groups[i].typ, typ)
		}
	}
	if len(groups[1].records) != 2 || groups[1].records[0] != start1 || groups[1].records[1] != start2 {
		t.Fatalf("start group = %v, want [start1,start2]", groups[1].records)
	}
	if len(groups[2].records) != 1 || groups[2].records[0] != finish {
		t.Fatalf("finish group = %v, want [finish]", groups[2].records)
	}
}

func TestGroupRecordsEmpty(t *testing.T) {
	if got := groupRecords(nil); got != nil {
		t.Fatalf("groupRecords(nil) = %v, want nil", got)
	}
}
