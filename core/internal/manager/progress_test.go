package manager

import (
	"testing"

	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func progressSave(num int64, policy savev1.SavePolicy, typ savev1.SaveType) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{Policy: policy, Type: typ},
		},
	}
}

func TestProgressTrackerContiguousWatermark(t *testing.T) {
	p := newProgressTracker(0)
	if got := p.watermark(); got != 0 {
		t.Fatalf("initial watermark = %d, want 0", got)
	}

	r1 := backfillScalar(1)
	r2 := backfillScalar(2)
	r3 := backfillScalar(3)
	p.onStored(r1, 1)
	p.onStored(r2, 2)
	p.onStored(r3, 3)

	// 先确认 seq=3，水位不能越过缺口 seq=1,2。
	p.resolve([]*recordv1.Record{r3})
	if got := p.watermark(); got != 0 {
		t.Fatalf("watermark with gap at 1 = %d, want 0", got)
	}

	// seq=1 被同 num 覆盖丢弃，此时 seq=2 仍未确认，水位推进到 1。
	p.resolveOne(r1)
	if got := p.watermark(); got != 1 {
		t.Fatalf("watermark after drop 1 = %d, want 1", got)
	}

	// seq=2 确认后，连续推进过 seq=2 和先前乱序确认的 seq=3。
	p.resolve([]*recordv1.Record{r2})
	if got := p.watermark(); got != 3 {
		t.Fatalf("watermark after gap fill = %d, want 3", got)
	}
}

func TestProgressTrackerNonUploadableAndUnknown(t *testing.T) {
	p := newProgressTracker(10)
	if got := p.watermark(); got != 10 {
		t.Fatalf("base watermark = %d, want 10", got)
	}

	// Start / END+CUSTOM / Finish 立即完结；CONFIG / NOW SAVE 需要上传确认。
	p.onStored(backfillStart(), 11)
	if got := p.watermark(); got != 11 {
		t.Fatalf("watermark after start = %d, want 11", got)
	}

	endSave := backfillEndSave(1)
	configSave := progressSave(2, savev1.SavePolicy_SAVE_POLICY_NOW, savev1.SaveType_SAVE_TYPE_CONFIG)
	s1 := backfillScalar(3)
	fin := backfillFinish()

	p.onStored(endSave, 12)
	p.onStored(configSave, 13)
	p.onStored(s1, 14)
	p.onStored(fin, 15)

	// 12(endSave) 自动完结但 13(configSave) 未确认，水位停在 12。
	if got := p.watermark(); got != 12 {
		t.Fatalf("watermark before configSave confirm = %d, want 12", got)
	}

	// 未知记录（例如仅上传不落盘的 crash error_log）不会影响任何状态。
	p.resolve([]*recordv1.Record{backfillScalar(999)})
	p.resolveOne(backfillScalar(998))
	if got := p.watermark(); got != 12 {
		t.Fatalf("watermark after unknown resolve = %d, want 12", got)
	}

	// 确认 14(s1) 与 13(configSave) 后，连续越过 15(fin)。
	p.resolve([]*recordv1.Record{s1, configSave})
	if got := p.watermark(); got != 15 {
		t.Fatalf("watermark after confirming 13,14 = %d, want 15", got)
	}
}

func TestProgressTrackerReplayedPointerReplacement(t *testing.T) {
	p := newProgressTracker(0)
	orig := backfillScalar(1)
	replayed := backfillScalar(1)

	p.onStored(orig, 1)
	// backfill 从本地文件重新反序列化出新指针替换旧指针。
	p.onStored(replayed, 1)

	// 旧指针即使被传入也已失效；只有重放后的新指针确认才能推进水位。
	p.resolve([]*recordv1.Record{orig})
	if got := p.watermark(); got != 0 {
		t.Fatalf("watermark after stale pointer = %d, want 0", got)
	}

	p.resolve([]*recordv1.Record{replayed})
	if got := p.watermark(); got != 1 {
		t.Fatalf("watermark after replayed pointer = %d, want 1", got)
	}
}

func TestProgressTrackerOutOfOrderLiveAndBackfill(t *testing.T) {
	// 模拟 resume 场景：base=100，live 先注册并确认 251，而 101..102 尚未由 backfill 扫描注册。
	p := newProgressTracker(100)
	live251 := backfillScalar(200)
	p.onStored(live251, 251)
	p.resolve([]*recordv1.Record{live251})

	// 101..250 尚未完结，水位决不能跳到 251。
	if got := p.watermark(); got != 100 {
		t.Fatalf("watermark jumped past unregistered gap: got %d, want 100", got)
	}

	r101 := backfillScalar(101)
	r102 := backfillScalar(102)
	p.onStored(r101, 101)
	p.onStored(r102, 102)
	p.resolve([]*recordv1.Record{r101, r102})
	if got := p.watermark(); got != 102 {
		t.Fatalf("watermark = %d, want 102", got)
	}
}
