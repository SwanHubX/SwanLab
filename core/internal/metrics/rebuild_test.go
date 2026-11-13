package metrics

import (
	"errors"
	"math"
	"os"
	"path/filepath"
	"testing"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func writeStoreRecords(t *testing.T, recs ...*recordv1.Record) (string, []int64) {
	t.Helper()
	path := filepath.Join(t.TempDir(), "run-test.swanlab")
	w, err := store.Create(path)
	if err != nil {
		t.Fatalf("Create: %v", err)
	}
	ends := make([]int64, len(recs))
	for i, rec := range recs {
		payload, merr := proto.Marshal(rec)
		if merr != nil {
			_ = w.Close()
			t.Fatalf("Marshal[%d]: %v", i, merr)
		}
		if werr := w.Write(payload); werr != nil {
			_ = w.Close()
			t.Fatalf("Write[%d]: %v", i, werr)
		}
		ends[i] = w.Offset()
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return path, ends
}

func colRecord(num int64, col *columnv1.ColumnRecord) *recordv1.Record {
	return &recordv1.Record{Num: num, RecordType: &recordv1.Record_Column{Column: col}}
}

func scalarDataRecord(num int64, key string, step int64, val float64) *recordv1.Record {
	return &recordv1.Record{Num: num, RecordType: &recordv1.Record_Scalar{Scalar: scalarRec(key, step, val)}}
}

func mediaDataRecord(num int64, key string, step int64) *recordv1.Record {
	return &recordv1.Record{Num: num, RecordType: &recordv1.Record_Media{Media: mediaRec(key, step)}}
}

func histogramDataRecord(num int64, key string, step int64) *recordv1.Record {
	return &recordv1.Record{Num: num, RecordType: &recordv1.Record_Histogram{Histogram: histogramRec(key, step)}}
}

func logDataRecord(num, epoch int64) *recordv1.Record {
	return &recordv1.Record{Num: num, RecordType: &recordv1.Record_Log{Log: &terminalv1.LogRecord{
		Line: "hello", Epoch: epoch,
	}}}
}

func saveDataRecord(num int64, name string, typ savev1.SaveType, policy savev1.SavePolicy) *recordv1.Record {
	return &recordv1.Record{Num: num, RecordType: &recordv1.Record_Save{Save: &savev1.SaveRecord{
		Name: name, Type: typ, Policy: policy,
	}}}
}

func TestRebuildEmptyHeaderOnly(t *testing.T) {
	path, _ := writeStoreRecords(t)
	res, err := Rebuild(path)
	if err != nil {
		t.Fatalf("Rebuild: %v", err)
	}
	if res.Records != 0 || res.EndOffset != store.HeaderLen ||
		res.Counter != 0 || res.Epoch != 0 ||
		res.GlobalStep != -1 || res.GlobalSystemStep != -1 ||
		res.Finished || res.HasStart || len(res.EndSaves) != 0 ||
		res.Metrics.Len() != 0 {
		t.Fatalf("unexpected empty rebuild result: %+v", res)
	}
}

func TestRebuildFullState(t *testing.T) {
	lossCol := builder.BuildAutoColumn(-1, "train/loss", columnv1.ColumnType_COLUMN_TYPE_SCALAR)
	sysCol := &columnv1.ColumnRecord{
		ColumnKey:   "__swanlab__.cpu",
		ColumnType:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		ColumnClass: columnv1.ColumnClass_COLUMN_CLASS_SYSTEM,
	}
	imgCol := builder.BuildAutoColumn(-1, "samples/img", columnv1.ColumnType_COLUMN_TYPE_IMAGE)
	histCol := builder.BuildAutoColumn(-1, "weights", columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM)

	recs := []*recordv1.Record{
		{Num: record.NumStart, RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}}},
		// 内部 config 保存使用负数槽位且即使配 END 也不应进 EndSaves。
		saveDataRecord(record.NumConfigSave, "config.yaml", savev1.SaveType_SAVE_TYPE_CONFIG, savev1.SavePolicy_SAVE_POLICY_END),
		colRecord(1, lossCol),
		colRecord(2, lossCol), // 重复列应静默跳过
		scalarDataRecord(3, "train/loss", 0, 2.5),
		scalarDataRecord(4, "train/loss", 0, 9.9), // 重复 step 应跳过，不更新极值
		scalarDataRecord(5, "train/loss", 1, math.NaN()),
		scalarDataRecord(6, "train/loss", 2, 1.25),
		colRecord(7, sysCol),
		scalarDataRecord(8, "__swanlab__.cpu", 10, 42.0),
		colRecord(9, imgCol),
		mediaDataRecord(10, "samples/img", 5),
		colRecord(11, histCol),
		histogramDataRecord(12, "weights", 7),
		logDataRecord(13, 3),
		logDataRecord(14, 9),
		saveDataRecord(15, "now.txt", savev1.SaveType_SAVE_TYPE_CUSTOM, savev1.SavePolicy_SAVE_POLICY_NOW),
		saveDataRecord(16, "end-1.pt", savev1.SaveType_SAVE_TYPE_CUSTOM, savev1.SavePolicy_SAVE_POLICY_END),
		saveDataRecord(18, "end-2.pt", savev1.SaveType_SAVE_TYPE_CUSTOM, savev1.SavePolicy_SAVE_POLICY_END),
		{Num: record.NumFinish, RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{
			State: runv1.RunState_RUN_STATE_FINISHED,
		}}},
	}

	path, ends := writeStoreRecords(t, recs...)
	res, err := Rebuild(path)
	if err != nil {
		t.Fatalf("Rebuild: %v", err)
	}

	if res.Records != uint64(len(recs)) || res.EndOffset != ends[len(ends)-1] {
		t.Fatalf("records/offset = (%d, %d), want (%d, %d)",
			res.Records, res.EndOffset, len(recs), ends[len(ends)-1])
	}
	if !res.HasStart || !res.Finished {
		t.Fatalf("HasStart=%v Finished=%v, want true/true", res.HasStart, res.Finished)
	}
	if res.Counter != 18 {
		t.Fatalf("Counter = %d, want 18", res.Counter)
	}
	if res.Epoch != 9 {
		t.Fatalf("Epoch = %d, want 9", res.Epoch)
	}
	if res.GlobalStep != 7 {
		t.Fatalf("GlobalStep = %d, want 7", res.GlobalStep)
	}
	if res.GlobalSystemStep != 10 {
		t.Fatalf("GlobalSystemStep = %d, want 10", res.GlobalSystemStep)
	}

	if len(res.EndSaves) != 2 ||
		res.EndSaves[0].GetSave().GetName() != "end-1.pt" ||
		res.EndSaves[1].GetSave().GetName() != "end-2.pt" {
		t.Fatalf("EndSaves = %+v, want [end-1.pt, end-2.pt]", res.EndSaves)
	}

	// 新 step 接纳与重复 step 去重验证。
	if p := res.Metrics.PrecheckScalar(scalarRec("train/loss", 2, 0.5)); p.NeedsColumn || p.Decision != DecisionSkip {
		t.Fatalf("duplicate step precheck = %+v, want NeedsColumn=false DecisionSkip", p)
	}
	if p := res.Metrics.PrecheckScalar(scalarRec("train/loss", 3, 0.5)); p.NeedsColumn || p.Decision != DecisionAccept {
		t.Fatalf("new step precheck = %+v, want NeedsColumn=false DecisionAccept", p)
	}
	if card := res.Metrics.cardinality("train/loss"); card != 3 {
		t.Fatalf("train/loss step cardinality = %d, want 3", card)
	}
	st := res.Metrics.metrics["train/loss"].scalar
	if st == nil || !st.has || st.min != 1.25 || st.max != 2.5 || st.latest != 1.25 {
		t.Fatalf("scalar extremes = %+v, want min=1.25 max=2.5 latest=1.25", st)
	}
}

func TestRebuildMissingColumnFallback(t *testing.T) {
	recs := []*recordv1.Record{
		scalarDataRecord(1, "acc", 4, 0.8),
		scalarDataRecord(2, "__swanlab__.mem", 12, 512),
		mediaDataRecord(3, "img", 6),
		histogramDataRecord(4, "grad", 9),
	}
	path, _ := writeStoreRecords(t, recs...)
	res, err := Rebuild(path)
	if err != nil {
		t.Fatalf("Rebuild: %v", err)
	}

	if res.Metrics.Len() != 4 {
		t.Fatalf("Metrics.Len() = %d, want 4", res.Metrics.Len())
	}
	if res.GlobalStep != 9 || res.GlobalSystemStep != 12 {
		t.Fatalf("GlobalStep=%d GlobalSystemStep=%d, want (9, 12)", res.GlobalStep, res.GlobalSystemStep)
	}

	// 验证降级列类型/类别与 resume 路径一致。
	if col := res.Metrics.metrics["acc"].column; col.ColumnType != columnv1.ColumnType_COLUMN_TYPE_SCALAR ||
		col.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_CUSTOM {
		t.Fatalf("fallback custom scalar column = %+v", col)
	}
	if col := res.Metrics.metrics["__swanlab__.mem"].column; col.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_SYSTEM {
		t.Fatalf("fallback system scalar column = %+v", col)
	}
	if col := res.Metrics.metrics["img"].column; col.ColumnType != columnv1.ColumnType_COLUMN_TYPE_UNSPECIFIED {
		t.Fatalf("fallback media column = %+v, want UNSPECIFIED", col)
	}
	if col := res.Metrics.metrics["grad"].column; col.ColumnType != columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM {
		t.Fatalf("fallback histogram column = %+v", col)
	}

	// 降级列随后可继续接收新数据。
	if ok, cerr := res.Metrics.CommitMedia(mediaRec("img", 7), nil); !ok || cerr != nil {
		t.Fatalf("CommitMedia after fallback = (%v, %v)", ok, cerr)
	}
	if ok, cerr := res.Metrics.CommitHistogram(histogramRec("grad", 10), nil); !ok || cerr != nil {
		t.Fatalf("CommitHistogram after fallback = (%v, %v)", ok, cerr)
	}
}

func TestRebuildStopsAtCorruptAndTruncatedTail(t *testing.T) {
	baseRecs := []*recordv1.Record{
		colRecord(1, builder.BuildAutoColumn(-1, "loss", columnv1.ColumnType_COLUMN_TYPE_SCALAR)),
		scalarDataRecord(2, "loss", 0, 1.0),
		scalarDataRecord(3, "loss", 1, 2.0),
	}

	t.Run("corrupt crc keeps valid prefix", func(t *testing.T) {
		path, ends := writeStoreRecords(t, baseRecs...)
		raw, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("ReadFile: %v", err)
		}
		// 破坏第 3 条记录末尾字节使 CRC 校验失败。
		raw[len(raw)-1] ^= 0xFF
		if werr := os.WriteFile(path, raw, 0o644); werr != nil {
			t.Fatalf("WriteFile: %v", werr)
		}

		res, rerr := Rebuild(path)
		if rerr != nil {
			t.Fatalf("Rebuild: %v", rerr)
		}
		if res.Records != 2 || res.EndOffset != ends[1] || res.GlobalStep != 0 || res.Counter != 2 {
			t.Fatalf("corrupt rebuild = %+v, want Records=2 EndOffset=%d GlobalStep=0 Counter=2", res, ends[1])
		}
	})

	t.Run("truncated tail keeps valid prefix and allows append", func(t *testing.T) {
		path, ends := writeStoreRecords(t, baseRecs...)
		raw, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("ReadFile: %v", err)
		}
		// 在第 3 条记录中间截断。
		if werr := os.WriteFile(path, raw[:ends[1]+5], 0o644); werr != nil {
			t.Fatalf("WriteFile: %v", werr)
		}

		res, rerr := Rebuild(path)
		if rerr != nil {
			t.Fatalf("Rebuild: %v", rerr)
		}
		if res.Records != 2 || res.EndOffset != ends[1] || res.GlobalStep != 0 || res.Counter != 2 {
			t.Fatalf("truncated rebuild = %+v, want Records=2 EndOffset=%d GlobalStep=0 Counter=2", res, ends[1])
		}

		// 验证 EndOffset 可直接传给 store.OpenAppend 追加新记录。
		w, oerr := store.OpenAppend(path, res.EndOffset)
		if oerr != nil {
			t.Fatalf("OpenAppend: %v", oerr)
		}
		p, merr := proto.Marshal(scalarDataRecord(3, "loss", 5, 0.5))
		if merr != nil {
			_ = w.Close()
			t.Fatalf("Marshal: %v", merr)
		}
		if werr := w.Write(p); werr != nil {
			_ = w.Close()
			t.Fatalf("Write: %v", werr)
		}
		if cerr := w.Close(); cerr != nil {
			t.Fatalf("Close: %v", cerr)
		}

		res2, rerr2 := Rebuild(path)
		if rerr2 != nil {
			t.Fatalf("Rebuild after append: %v", rerr2)
		}
		if res2.Records != 3 || res2.GlobalStep != 5 || res2.Counter != 3 {
			t.Fatalf("rebuild after append = %+v, want Records=3 GlobalStep=5 Counter=3", res2)
		}
	})

	t.Run("corrupt protobuf payload stops safely", func(t *testing.T) {
		path, ends := writeStoreRecords(t, baseRecs[:2]...)
		w, oerr := store.OpenAppend(path, ends[1])
		if oerr != nil {
			t.Fatalf("OpenAppend: %v", oerr)
		}
		// 写入合法 chunk 帧但非法的 protobuf wire 字节。
		if werr := w.Write([]byte{0xFF, 0xFF, 0xFF, 0xFF}); werr != nil {
			_ = w.Close()
			t.Fatalf("Write: %v", werr)
		}
		if cerr := w.Close(); cerr != nil {
			t.Fatalf("Close: %v", cerr)
		}

		res, rerr := Rebuild(path)
		if rerr != nil {
			t.Fatalf("Rebuild: %v", rerr)
		}
		if res.Records != 2 || res.EndOffset != ends[1] {
			t.Fatalf("corrupt proto rebuild = %+v, want Records=2 EndOffset=%d", res, ends[1])
		}
	})
}

func TestRebuildOpenError(t *testing.T) {
	if _, err := Rebuild(filepath.Join(t.TempDir(), "missing.swanlab")); !errors.Is(err, os.ErrNotExist) {
		t.Fatalf("missing file err = %v, want ErrNotExist", err)
	}
	bad := filepath.Join(t.TempDir(), "bad.swanlab")
	if err := os.WriteFile(bad, []byte("bad"), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	if _, err := Rebuild(bad); !errors.Is(err, store.ErrBadHeader) {
		t.Fatalf("bad header err = %v, want ErrBadHeader", err)
	}
}
