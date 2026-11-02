package metrics

import (
	"errors"
	"math"
	"testing"

	"github.com/RoaringBitmap/roaring/v2/roaring64"

	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
)

func scalarRec(key string, step int64, value float64) *datav1.ScalarRecord {
	return &datav1.ScalarRecord{
		Key: key, Step: step, Type: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		Value: &datav1.ScalarValue{Number: value},
	}
}

func mediaRec(key string, step int64) *datav1.MediaRecord {
	return &datav1.MediaRecord{Key: key, Step: step, Type: columnv1.ColumnType_COLUMN_TYPE_IMAGE}
}

func autoCol(key string) *columnv1.ColumnRecord {
	return &columnv1.ColumnRecord{
		ColumnKey:  key,
		ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
	}
}

// cardinality 断言用的位图计数（白盒）。
func (m *RunMetrics) cardinality(key string) uint64 {
	return m.metrics[key].steps.GetCardinality()
}

func TestDefineDuplicate(t *testing.T) {
	m := New()
	if err := m.DefineScalar("loss", autoCol("loss"), -1); err != nil {
		t.Fatal(err)
	}
	err := m.DefineScalar("loss", autoCol("loss"), -1)
	if !errors.Is(err, ErrColumnExists) {
		t.Fatalf("重复定义应报 ErrColumnExists, got %v", err)
	}
	if err := m.DefineMedia("pic", nil, -1); !errors.Is(err, ErrNilColumn) {
		t.Fatalf("nil 列定义应报 ErrNilColumn, got %v", err)
	}
	if err := m.DefineMedia("", autoCol(""), -1); !errors.Is(err, ErrEmptyKey) {
		t.Fatalf("空 key 应报 ErrEmptyKey, got %v", err)
	}
	if !m.Defined("loss") || m.Defined("pic") {
		t.Fatal("Defined 状态错误")
	}
	if m.Len() != 1 {
		t.Fatalf("Len = %d, want 1", m.Len())
	}
	if ms, ok := m.MinStep("loss"); !ok || ms != -1 {
		t.Fatalf("MinStep = (%d,%v), want (-1,true)", ms, ok)
	}
	if ms, ok := m.MinStep("nope"); ok || ms != 0 {
		t.Fatalf("未定义 key 的 MinStep = (%d,%v), want (0,false)", ms, ok)
	}
}

func TestPrecheckIsReadOnly(t *testing.T) {
	// 预检失败不写位图：用计数断言（完成定义的核心用例）
	m := New()
	if err := m.DefineScalar("loss", autoCol("loss"), -1); err != nil {
		t.Fatal(err)
	}

	cases := []struct {
		name string
		rec  *datav1.ScalarRecord
	}{
		{"重复 step", scalarRec("loss", 1, 0.1)},
		{"过期 step", scalarRec("loss", -1, 0.1)},
		{"超限 step", scalarRec("loss", MaxStep+1, 0.1)},
		{"类型冲突", &datav1.ScalarRecord{Key: "loss", Step: 1, Type: columnv1.ColumnType_COLUMN_TYPE_IMAGE}},
		{"新 key 超限", scalarRec("new", MaxStep*2, 0.1)},
		{"新 key 过期", scalarRec("new", -5, 0.1)},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			before := m.cardinalitySnapshot()
			p := m.PrecheckScalar(tc.rec)
			after := m.cardinalitySnapshot()
			if before != after {
				t.Fatalf("预检修改了位图: %v -> %v", before, after)
			}
			if p.Err == nil && p.Decision == DecisionAccept {
				// 合法接纳的预检也不允许立即写位图
				if m.cardinality("loss") != 0 {
					t.Fatal("预检阶段写入了位图")
				}
			}
		})
	}
}

// cardinalitySnapshot 汇总全部指标的位图计数，检测任何写入。
func (m *RunMetrics) cardinalitySnapshot() uint64 {
	var total uint64
	for _, mt := range m.metrics {
		total += mt.steps.GetCardinality()
	}
	return total
}

func TestPrecheckFreshKeyNeedsColumn(t *testing.T) {
	m := New()
	// 新 key 且 step 合法：NeedsColumn + Accept
	p := m.PrecheckScalar(scalarRec("loss", 0, 0.1))
	if !p.NeedsColumn || p.Decision != DecisionAccept || p.Err != nil {
		t.Fatalf("Precheck = %+v, want NeedsColumn+Accept", p)
	}
	// 新 key 但 step 超限：列仍会定义，数据拒绝。
	p = m.PrecheckScalar(scalarRec("loss", MaxStep+1, 0.1))
	if !p.NeedsColumn || p.Decision != DecisionSkip || !errors.Is(p.Err, ErrStepLimit) {
		t.Fatalf("Precheck = %+v, want NeedsColumn+Skip+ErrStepLimit", p)
	}
	// 新 key 且 step ≤ -1：列定义 + 数据跳过
	p = m.PrecheckScalar(scalarRec("loss", -1, 0.1))
	if !p.NeedsColumn || p.Decision != DecisionSkip || p.Err != nil {
		t.Fatalf("Precheck = %+v, want NeedsColumn+Skip", p)
	}
}

func TestCommitScalarLifecycle(t *testing.T) {
	m := New()
	c := autoCol("loss")

	// step 0 是合法首步（min_step=-1 ⇒ step ≤ -1 拒绝）
	if ok, err := m.CommitScalar(scalarRec("loss", 0, 1.5), c); !ok || err != nil {
		t.Fatalf("首次提交 = (%v, %v)", ok, err)
	}
	if m.cardinality("loss") != 1 {
		t.Fatalf("位图计数 = %d, want 1", m.cardinality("loss"))
	}
	// 重复 step：静默跳过，不报错
	if ok, err := m.CommitScalar(scalarRec("loss", 0, 9.9), nil); ok || err != nil {
		t.Fatalf("重复 step = (%v, %v), want (false, nil)", ok, err)
	}
	// 负 step：水位线拒绝
	if ok, err := m.CommitScalar(scalarRec("loss", -1, 9.9), nil); ok || err != nil {
		t.Fatalf("负 step = (%v, %v), want (false, nil)", ok, err)
	}
	// 超限 step：ErrStepLimit
	if _, err := m.CommitScalar(scalarRec("loss", MaxStep+1, 1.0), nil); !errors.Is(err, ErrStepLimit) {
		t.Fatalf("超限 step err = %v, want ErrStepLimit", err)
	}
	// 类型冲突：ErrTypeMismatch
	if _, err := m.CommitScalar(&datav1.ScalarRecord{
		Key: "loss", Step: 1,
		Type: columnv1.ColumnType_COLUMN_TYPE_IMAGE,
	}, nil); !errors.Is(err, ErrTypeMismatch) {
		t.Fatalf("类型冲突 err = %v, want ErrTypeMismatch", err)
	}
	// 位图不受失败提交影响
	if m.cardinality("loss") != 1 {
		t.Fatalf("失败提交污染位图: %d", m.cardinality("loss"))
	}
	// 流程错误：已定义 key 携带列定义 / 未定义 key 缺列定义
	if _, err := m.CommitScalar(scalarRec("loss", 2, 1.0), autoCol("loss")); !errors.Is(err, ErrColumnExists) {
		t.Fatalf("已定义携带列 = %v, want ErrColumnExists", err)
	}
	if _, err := m.CommitScalar(scalarRec("ghost", 2, 1.0), nil); !errors.Is(err, ErrColumnMissing) {
		t.Fatalf("新 key 缺列 = %v, want ErrColumnMissing", err)
	}
	if m.Defined("ghost") {
		t.Fatal("ErrColumnMissing 不应定义新列")
	}
}

func TestStepBoundaryOrder(t *testing.T) {
	// 判定顺序：min_step 先于超限判定。
	m := New()
	if err := m.DefineScalar("k", autoCol("k"), MaxStep+10); err != nil {
		t.Fatal(err)
	}
	if _, err := m.CommitScalar(scalarRec("k", MaxStep+5, 1.0), nil); err != nil {
		t.Fatalf("水位线之上的病态 step 应按 Skip 处理, got %v", err)
	}
	// 恰好等于水位线的 step 拒绝
	if ok, _ := m.CommitScalar(scalarRec("k", MaxStep+10, 1.0), nil); ok {
		t.Fatal("step == min_step 应拒绝")
	}

	// 常规水位线：== 拒绝、+1 接纳
	m2 := New()
	if err := m2.DefineScalar("k", autoCol("k"), 49); err != nil {
		t.Fatal(err)
	}
	if ok, _ := m2.CommitScalar(scalarRec("k", 49, 1.0), nil); ok {
		t.Fatal("step == min_step 应拒绝")
	}
	if ok, err := m2.CommitScalar(scalarRec("k", 50, 1.0), nil); !ok || err != nil {
		t.Fatalf("step = min_step+1 = (%v, %v), want 接纳", ok, err)
	}

	// 恰好等于 MaxStep 的 step 合法
	m3 := New()
	if ok, err := m3.CommitScalar(scalarRec("k", MaxStep, 1.0), autoCol("k")); !ok || err != nil {
		t.Fatalf("step == MaxStep = (%v, %v)", ok, err)
	}
}

func TestScalarExtremesFiniteOnly(t *testing.T) {
	m := New()
	c := autoCol("loss")

	// 有限值推进极值
	if ok, err := m.CommitScalar(scalarRec("loss", 0, 1.0), c); !ok || err != nil {
		t.Fatal(err)
	}
	if ok, err := m.CommitScalar(scalarRec("loss", 1, -2.0), nil); !ok || err != nil {
		t.Fatalf("commit loss step 1 = (%v, %v)", ok, err)
	}
	if ok, err := m.CommitScalar(scalarRec("loss", 2, 3.0), nil); !ok || err != nil {
		t.Fatalf("commit loss step 2 = (%v, %v)", ok, err)
	}
	st := m.metrics["loss"].scalar
	if !st.has || st.latest != 3.0 || st.min != -2.0 || st.max != 3.0 {
		t.Fatalf("极值 = latest %v min %v max %v has %v", st.latest, st.min, st.max, st.has)
	}
	// NaN / ±Inf 不更新极值
	for i, v := range []float64{math.NaN(), math.Inf(1), math.Inf(-1)} {
		if ok, err := m.CommitScalar(scalarRec("loss", int64(3+i), v), nil); !ok || err != nil {
			t.Fatalf("commit loss step %d = (%v, %v)", 3+i, ok, err)
		}
		if st.latest != 3.0 || st.min != -2.0 || st.max != 3.0 {
			t.Fatalf("非法值 %v 更新了极值: %+v", v, st)
		}
	}
	// NaN/Inf 的 step 仍登记位图（step 与值合法性解耦）
	if m.cardinality("loss") != 6 {
		t.Fatalf("位图计数 = %d, want 6", m.cardinality("loss"))
	}
}

// TestCommitNewKeyRejectedStepStillDefinesColumn：新 key 因 step 非法被拒绝时列仍已定义，
// 调用方必须继续投递该列记录（对齐 Python safe.block）。
func TestCommitNewKeyRejectedStepStillDefinesColumn(t *testing.T) {
	m := New()
	col := autoCol("hot")
	if ok, err := m.CommitScalar(scalarRec("hot", MaxStep+1, 1.0), col); ok || !errors.Is(err, ErrStepLimit) {
		t.Fatalf("新 key 超限 = (%v, %v), want (false, ErrStepLimit)", ok, err)
	}
	if !m.Defined("hot") {
		t.Fatal("被拒绝的超限记录也应定义新列")
	}
	if ms, ok := m.MinStep("hot"); !ok || ms != -1 {
		t.Fatalf("MinStep = (%d,%v), want (-1,true)", ms, ok)
	}
	// 列已定义后，合法 step 不再需要列记录
	if ok, err := m.CommitScalar(scalarRec("hot", 0, 1.0), nil); !ok || err != nil {
		t.Fatalf("合法 step = (%v, %v), want 接纳", ok, err)
	}
	if m.cardinality("hot") != 1 {
		t.Fatalf("位图计数 = %d, want 1", m.cardinality("hot"))
	}
}

func TestCommitMediaAndTypeMatch(t *testing.T) {
	m := New()
	// 新媒体列 + IMAGE 类型
	if ok, err := m.CommitMedia(mediaRec("pic", 0), &columnv1.ColumnRecord{
		ColumnKey: "pic", ColumnType: columnv1.ColumnType_COLUMN_TYPE_IMAGE,
	}); !ok || err != nil {
		t.Fatalf("首次媒体提交 = (%v, %v)", ok, err)
	}
	// 已定义 IMAGE 列收到 AUDIO：类型冲突
	if _, err := m.CommitMedia(&datav1.MediaRecord{
		Key: "pic", Step: 1,
		Type: columnv1.ColumnType_COLUMN_TYPE_AUDIO,
	}, nil); !errors.Is(err, ErrTypeMismatch) {
		t.Fatalf("媒体类型冲突 err = %v", err)
	}
	// 重复 step 跳过
	if ok, err := m.CommitMedia(mediaRec("pic", 0), nil); ok || err != nil {
		t.Fatalf("重复媒体 step = (%v, %v)", ok, err)
	}
	// 标量落到 IMAGE 列：类型冲突
	if _, err := m.CommitScalar(scalarRec("pic", 2, 1.0), nil); !errors.Is(err, ErrTypeMismatch) {
		t.Fatalf("标量落媒体列 err = %v", err)
	}
}

func TestCommitMediaFlowErrors(t *testing.T) {
	m := New()
	// 新 key 未携带列定义：流程错误，且不得定义列
	if _, err := m.CommitMedia(mediaRec("ghost", 0), nil); !errors.Is(err, ErrColumnMissing) {
		t.Fatalf("新 key 缺列 = %v, want ErrColumnMissing", err)
	}
	if m.Defined("ghost") {
		t.Fatal("ErrColumnMissing 不应定义新列")
	}
	// 已定义 key 再次携带列定义：流程错误
	col := &columnv1.ColumnRecord{ColumnKey: "pic", ColumnType: columnv1.ColumnType_COLUMN_TYPE_IMAGE}
	if ok, err := m.CommitMedia(mediaRec("pic", 0), col); !ok || err != nil {
		t.Fatalf("首次提交 = (%v, %v)", ok, err)
	}
	if _, err := m.CommitMedia(mediaRec("pic", 1), col); !errors.Is(err, ErrColumnExists) {
		t.Fatalf("已定义携带列 = %v, want ErrColumnExists", err)
	}
}

func TestPrecheckMedia(t *testing.T) {
	m := New()
	// 新 key：媒体列由调用方构建，step 合法则接纳
	if p := m.PrecheckMedia(mediaRec("pic", 0)); !p.NeedsColumn || p.Decision != DecisionAccept || p.Err != nil {
		t.Fatalf("新 key 预检 = %+v, want NeedsColumn+Accept", p)
	}
	// 新 key 且 step 超限：列仍会定义，数据拒绝
	if p := m.PrecheckMedia(mediaRec("pic", MaxStep+1)); !p.NeedsColumn ||
		p.Decision != DecisionSkip || !errors.Is(p.Err, ErrStepLimit) {
		t.Fatalf("新 key 超限预检 = %+v, want NeedsColumn+Skip+ErrStepLimit", p)
	}
	if m.Len() != 0 {
		t.Fatal("预检不应定义列")
	}

	col := &columnv1.ColumnRecord{ColumnKey: "pic", ColumnType: columnv1.ColumnType_COLUMN_TYPE_IMAGE}
	if ok, err := m.CommitMedia(mediaRec("pic", 0), col); !ok || err != nil {
		t.Fatalf("提交 = (%v, %v)", ok, err)
	}
	// 已定义列：重复 step 静默跳过，新 step 接纳
	if p := m.PrecheckMedia(mediaRec("pic", 0)); p.NeedsColumn || p.Decision != DecisionSkip || p.Err != nil {
		t.Fatalf("重复 step 预检 = %+v, want Skip", p)
	}
	if p := m.PrecheckMedia(mediaRec("pic", 1)); p.NeedsColumn || p.Decision != DecisionAccept || p.Err != nil {
		t.Fatalf("新 step 预检 = %+v, want Accept", p)
	}
	// 类型冲突：确定性错误，且不写位图
	audio := &datav1.MediaRecord{Key: "pic", Step: 1, Type: columnv1.ColumnType_COLUMN_TYPE_AUDIO}
	if p := m.PrecheckMedia(audio); p.Decision != DecisionSkip || !errors.Is(p.Err, ErrTypeMismatch) {
		t.Fatalf("类型冲突预检 = %+v, want Skip+ErrTypeMismatch", p)
	}
	if m.cardinality("pic") != 1 {
		t.Fatalf("预检改动了位图: %d, want 1", m.cardinality("pic"))
	}
}

func TestPrecheckMediaResumePlaceholder(t *testing.T) {
	// UNSPECIFIED 占位列放行任意媒体类型，但水位线依旧生效
	res, err := NewFromSummary(&Summary{Media: []SummaryItem{{Key: "camera", Step: 9}}})
	if err != nil {
		t.Fatal(err)
	}
	m := res.Metrics
	if p := m.PrecheckMedia(&datav1.MediaRecord{
		Key: "camera", Step: 9, Type: columnv1.ColumnType_COLUMN_TYPE_VIDEO,
	}); p.Decision != DecisionSkip || p.Err != nil {
		t.Fatalf("水位线预检 = %+v, want Skip", p)
	}
	if p := m.PrecheckMedia(&datav1.MediaRecord{
		Key: "camera", Step: 10, Type: columnv1.ColumnType_COLUMN_TYPE_VIDEO,
	}); p.NeedsColumn || p.Decision != DecisionAccept || p.Err != nil {
		t.Fatalf("占位列预检 = %+v, want Accept", p)
	}
}

func TestBitmapOptimizeCadence(t *testing.T) {
	// 写入 4096 个连续 step：越过压缩周期后位图语义不变，
	// 且 run 容器收敛（序列化体积远小于逐点数组）。
	m := New()
	if ok, err := m.CommitScalar(scalarRec("k", 0, 1.0), autoCol("k")); !ok || err != nil {
		t.Fatal(err)
	}
	for step := int64(1); step <= 4100; step++ {
		if ok, err := m.CommitScalar(scalarRec("k", step, 1.0), nil); !ok || err != nil {
			t.Fatalf("step %d = (%v, %v)", step, ok, err)
		}
	}
	bm := m.metrics["k"].steps
	if got := bm.GetCardinality(); got != 4101 {
		t.Fatalf("cardinality = %d, want 4101", got)
	}
	if !bm.Contains(0) || !bm.Contains(4096) || !bm.Contains(4100) {
		t.Fatal("位图成员缺失")
	}
	if sz := bm.GetSerializedSizeInBytes(); sz > 4096*2/8+1024 {
		t.Fatalf("run 容器未收敛, 序列化 %d 字节", sz)
	}
	// 乱序稀疏写入也安全（array 容器无收益时保持原状）
	m2 := New()
	if err := m2.DefineScalar("s", autoCol("s"), -1); err != nil {
		t.Fatal(err)
	}
	bm2 := roaring64.New()
	for i := int64(0); i < 100; i++ {
		step := i * 9973 % MaxStep
		if ok, err := m2.CommitScalar(scalarRec("s", step, 1.0), nil); !ok || err != nil {
			t.Fatalf("commit s step %d = (%v, %v)", step, ok, err)
		}
		bm2.Add(uint64(step))
	}
	if m2.metrics["s"].steps.GetCardinality() != bm2.GetCardinality() {
		t.Fatal("稀疏位图计数不一致")
	}
}

func TestTwoPhasePlanMatchesCommit(t *testing.T) {
	// 预检结论与提交结果成对出现：Accept ⇒ Commit 接纳；Skip/Err ⇒ Commit 对应拒绝
	m := New()
	plan := m.PrecheckScalar(scalarRec("a", 5, 0.5))
	var col *columnv1.ColumnRecord
	if plan.NeedsColumn {
		col = autoCol("a")
	}
	ok, err := m.CommitScalar(scalarRec("a", 5, 0.5), col)
	if plan.Decision == DecisionAccept && plan.Err == nil {
		if !ok || err != nil {
			t.Fatalf("预检 Accept 但提交 = (%v, %v)", ok, err)
		}
	} else if ok {
		t.Fatal("预检非 Accept 但提交接纳")
	}
}
