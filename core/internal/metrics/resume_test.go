package metrics

import (
	"errors"
	"testing"

	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
)

func TestNewFromSummaryNil(t *testing.T) {
	// 无摘要（新实验）：全默认值，等价于 New()
	res, err := NewFromSummary(nil)
	if err != nil {
		t.Fatal(err)
	}
	if res.ConsoleEpoch != 0 || res.GlobalStep != -1 || res.GlobalSystemStep != -1 {
		t.Fatalf("默认值 = (%d, %d, %d), want (0, -1, -1)",
			res.ConsoleEpoch, res.GlobalStep, res.GlobalSystemStep)
	}
	if res.Metrics.Len() != 0 {
		t.Fatal("空状态机应无指标")
	}

	res2, err := NewFromSummary(&Summary{})
	if err != nil {
		t.Fatal(err)
	}
	if res2.GlobalStep != -1 || res2.GlobalSystemStep != -1 || res2.ConsoleEpoch != 0 {
		t.Fatalf("空摘要默认值错误: %+v", res2)
	}
}

func TestNewFromSummaryFull(t *testing.T) {
	s := &Summary{
		Log:   []SummaryItem{{Key: "stdout", Step: 42}},
		Media: []SummaryItem{{Key: "pic", Step: 100}},
		Scalar: []SummaryItem{
			{Key: "loss", Step: 88},
			{Key: "__swanlab__.gpu.util", Step: 66},
			{Key: "acc", Step: 12},
		},
	}
	res, err := NewFromSummary(s)
	if err != nil {
		t.Fatal(err)
	}
	if res.ConsoleEpoch != 42 {
		t.Fatalf("ConsoleEpoch = %d, want 42", res.ConsoleEpoch)
	}
	// 媒体 100 与标量最大 88 竞争 GlobalStep；系统列只推进 GlobalSystemStep
	if res.GlobalStep != 100 {
		t.Fatalf("GlobalStep = %d, want 100", res.GlobalStep)
	}
	if res.GlobalSystemStep != 66 {
		t.Fatalf("GlobalSystemStep = %d, want 66", res.GlobalSystemStep)
	}

	m := res.Metrics
	if m.Len() != 4 {
		t.Fatalf("恢复指标数 = %d, want 4", m.Len())
	}
	// 水位线就位
	for key, want := range map[string]int64{"pic": 100, "loss": 88, "__swanlab__.gpu.util": 66, "acc": 12} {
		got, ok := m.MinStep(key)
		if !ok || got != want {
			t.Fatalf("MinStep(%q) = (%d,%v), want %d", key, got, ok, want)
		}
	}
	// 媒体恢复列是 UNSPECIFIED 占位，任意媒体类型放行
	if col := m.metrics["pic"].column; col.ColumnType != columnv1.ColumnType_COLUMN_TYPE_UNSPECIFIED {
		t.Fatalf("媒体恢复列类型 = %v, want UNSPECIFIED", col.ColumnType)
	}
	if _, err := m.CommitMedia(&datav1.MediaRecord{
		Key: "pic", Step: 101,
		Type: columnv1.ColumnType_COLUMN_TYPE_VIDEO,
	}, nil); err != nil {
		t.Fatalf("占位列收到 VIDEO 媒体应放行: %v", err)
	}
	// 标量恢复列类型 SCALAR；系统列类别 SYSTEM、自定义列类别 CUSTOM
	if col := m.metrics["loss"].column; col.ColumnType != columnv1.ColumnType_COLUMN_TYPE_SCALAR ||
		col.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_CUSTOM {
		t.Fatalf("自定义标量恢复列 = %+v", col)
	}
	if col := m.metrics["__swanlab__.gpu.util"].column; col.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_SYSTEM {
		t.Fatalf("系统恢复列类别 = %v", col.ColumnClass)
	}
}

func TestResumeWatermarkRejectsHistory(t *testing.T) {
	// 恢复后：≤ 水位线的 step 全部拒绝，水位线 +1 起接纳；位图从空开始不重建历史
	res, err := NewFromSummary(&Summary{Scalar: []SummaryItem{{Key: "loss", Step: 100}}})
	if err != nil {
		t.Fatal(err)
	}
	m := res.Metrics
	for _, step := range []int64{99, 100} {
		if ok, err := m.CommitScalar(scalarRec("loss", step, 1.0), nil); ok || err != nil {
			t.Fatalf("恢复后 step %d = (%v, %v), want 拒绝", step, ok, err)
		}
	}
	if ok, err := m.CommitScalar(scalarRec("loss", 101, 1.0), nil); !ok || err != nil {
		t.Fatalf("恢复后 step 101 = (%v, %v), want 接纳", ok, err)
	}
	if m.cardinality("loss") != 1 {
		t.Fatalf("恢复位图计数 = %d, want 1（仅新接纳的 step）", m.cardinality("loss"))
	}
	// 预检同样尊重水位线
	if p := m.PrecheckScalar(scalarRec("loss", 100, 1.0)); p.Decision != DecisionSkip || p.Err != nil {
		t.Fatalf("预检 = %+v, want Skip", p)
	}
}

func TestNewFromSummaryDuplicateKey(t *testing.T) {
	s := &Summary{Scalar: []SummaryItem{{Key: "loss", Step: 1}, {Key: "loss", Step: 2}}}
	if _, err := NewFromSummary(s); err == nil {
		t.Fatal("摘要重复 key 应报错")
	} else if !errors.Is(err, ErrColumnExists) {
		t.Fatalf("err = %v, want wrap ErrColumnExists", err)
	}
}

func TestNewFromSummaryNegativeSteps(t *testing.T) {
	// 云端摘要含负 step（异常数据）：水位线为负 ⇒ step 0 起可写
	res, err := NewFromSummary(&Summary{Scalar: []SummaryItem{{Key: "loss", Step: -3}}})
	if err != nil {
		t.Fatal(err)
	}
	if ok, err := res.Metrics.CommitScalar(scalarRec("loss", -2, 1.0), nil); !ok || err != nil {
		t.Fatalf("step -2 = (%v, %v), want 接纳（水位线 -3）", ok, err)
	}
	if ok, _ := res.Metrics.CommitScalar(scalarRec("loss", -3, 1.0), nil); ok {
		t.Fatal("step == 水位线应拒绝")
	}
}

func TestSystemKey(t *testing.T) {
	if !SystemKey("__swanlab__.gpu") || SystemKey("loss") || SystemKey("__swanlab__") {
		t.Fatal("SystemKey 前缀判定错误")
	}
}
