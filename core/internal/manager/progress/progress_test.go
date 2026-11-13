package progress

import (
	"testing"

	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func progressScalar(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Scalar{Scalar: &datav1.ScalarRecord{
			Key:   "loss",
			Step:  num,
			Type:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
			Value: &datav1.ScalarValue{Number: float64(num)},
		}},
	}
}

func progressSave(num int64, policy savev1.SavePolicy, typ savev1.SaveType) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{Policy: policy, Type: typ},
		},
	}
}

func TestProgressTrackerContiguousWatermark(t *testing.T) {
	p := NewTracker(0)
	if got := p.Watermark(); got != 0 {
		t.Fatalf("initial watermark = %d, want 0", got)
	}

	r1 := progressScalar(1)
	r2 := progressScalar(2)
	r3 := progressScalar(3)
	p.OnStored(r1, 1)
	p.OnStored(r2, 2)
	p.OnStored(r3, 3)

	// 先确认 seq=3，水位不能越过缺口 seq=1,2。
	p.Resolve([]*recordv1.Record{r3})
	if got := p.Watermark(); got != 0 {
		t.Fatalf("watermark with gap at 1 = %d, want 0", got)
	}

	// seq=1 被同 num 覆盖丢弃，此时 seq=2 仍未确认，水位推进到 1。
	p.ResolveOne(r1)
	if got := p.Watermark(); got != 1 {
		t.Fatalf("watermark after drop 1 = %d, want 1", got)
	}

	// seq=2 确认后，连续推进过 seq=2 和先前乱序确认的 seq=3。
	p.Resolve([]*recordv1.Record{r2})
	if got := p.Watermark(); got != 3 {
		t.Fatalf("watermark after gap fill = %d, want 3", got)
	}
}

func TestProgressTrackerNonUploadableAndUnknown(t *testing.T) {
	p := NewTracker(10)
	if got := p.Watermark(); got != 10 {
		t.Fatalf("base watermark = %d, want 10", got)
	}

	startRec := &recordv1.Record{Num: record.NumStart, RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}}}
	finishRec := &recordv1.Record{Num: record.NumFinish, RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{}}}

	// Start / END+CUSTOM / Finish 立即完结；CONFIG / NOW SAVE 需要上传确认。
	p.OnStored(startRec, 11)
	if got := p.Watermark(); got != 11 {
		t.Fatalf("watermark after start = %d, want 11", got)
	}

	endSave := progressSave(1, savev1.SavePolicy_SAVE_POLICY_END, savev1.SaveType_SAVE_TYPE_CUSTOM)
	configSave := progressSave(2, savev1.SavePolicy_SAVE_POLICY_NOW, savev1.SaveType_SAVE_TYPE_CONFIG)
	s1 := progressScalar(3)

	p.OnStored(endSave, 12)
	p.OnStored(configSave, 13)
	p.OnStored(s1, 14)
	p.OnStored(finishRec, 15)

	// 12(endSave) 自动完结但 13(configSave) 未确认，水位停在 12。
	if got := p.Watermark(); got != 12 {
		t.Fatalf("watermark before configSave confirm = %d, want 12", got)
	}

	// 未知记录（例如仅上传不落盘的 crash error_log）不会影响任何状态。
	p.Resolve([]*recordv1.Record{progressScalar(999)})
	p.ResolveOne(progressScalar(998))
	if got := p.Watermark(); got != 12 {
		t.Fatalf("watermark after unknown resolve = %d, want 12", got)
	}

	// 确认 14(s1) 与 13(configSave) 后，连续越过 15(fin)。
	p.Resolve([]*recordv1.Record{s1, configSave})
	if got := p.Watermark(); got != 15 {
		t.Fatalf("watermark after confirming 13,14 = %d, want 15", got)
	}
}

func TestProgressTrackerReplayedPointerReplacement(t *testing.T) {
	p := NewTracker(0)
	orig := progressScalar(1)
	replayed := progressScalar(1)

	p.OnStored(orig, 1)
	// backfill 从本地文件重新反序列化出新指针替换旧指针。
	p.OnStored(replayed, 1)

	// 旧指针即使被传入也已失效；只有重放后的新指针确认才能推进水位。
	p.Resolve([]*recordv1.Record{orig})
	if got := p.Watermark(); got != 0 {
		t.Fatalf("watermark after stale pointer = %d, want 0", got)
	}

	p.Resolve([]*recordv1.Record{replayed})
	if got := p.Watermark(); got != 1 {
		t.Fatalf("watermark after replayed pointer = %d, want 1", got)
	}
}

func TestProgressTrackerOutOfOrderLiveAndBackfill(t *testing.T) {
	// 模拟 resume 场景：base=100，live 先注册并确认 251，而 101..102 尚未由 backfill 扫描注册。
	p := NewTracker(100)
	live251 := progressScalar(200)
	p.OnStored(live251, 251)
	p.Resolve([]*recordv1.Record{live251})

	// 101..250 尚未完结，水位决不能跳到 251。
	if got := p.Watermark(); got != 100 {
		t.Fatalf("watermark jumped past unregistered gap: got %d, want 100", got)
	}

	r101 := progressScalar(101)
	r102 := progressScalar(102)
	p.OnStored(r101, 101)
	p.OnStored(r102, 102)
	p.Resolve([]*recordv1.Record{r101, r102})
	if got := p.Watermark(); got != 102 {
		t.Fatalf("watermark = %d, want 102", got)
	}
}
