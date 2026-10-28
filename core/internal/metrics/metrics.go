// Package metrics 实现单个 run 的指标状态机：标量与媒体列定义、基于
// roaring64 的 step 去重、有限浮点极值追踪，以及从云端摘要恢复状态。
//
// 状态机不加锁，由 Session 的写锁保护。接收分为两步，以配合上层的容量
// 预留与背压：
//
//	PrecheckScalar / PrecheckMedia  只读预检，绝不修改状态；
//	CommitScalar  / CommitMedia     正式提交，定义新列、登记 step、更新极值。
//
// 位图仅驻留内存，不持久化。
package metrics

import (
	"errors"
	"fmt"
	"math"

	"github.com/RoaringBitmap/roaring/v2/roaring64"

	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
)

// MaxStep 是 step 上限：超过此值视为病态输入，返回 ErrStepLimit。
const MaxStep int64 = 1 << 48

// 确定性错误哨兵。调用方（Session）把它们映射为「跳过该条记录」，
// 不影响批次内其余记录，也不参与任何重试。
var (
	// ErrStepLimit: step 超过 MaxStep。
	ErrStepLimit = errors.New("step exceeds the supported step limit")
	// ErrTypeMismatch: 已定义列与记录类型冲突。
	ErrTypeMismatch = errors.New("metric type mismatch")
	// ErrColumnExists: DefineXxx / Commit 携带列定义时 key 已存在。
	ErrColumnExists = errors.New("metric already defined")
	// ErrColumnMissing: 新 key 提交时未携带列定义（调用方流程错误）。
	ErrColumnMissing = errors.New("column record required for new metric key")
	// ErrEmptyKey / ErrNilColumn: DefineXxx 的入参约束。
	ErrEmptyKey  = errors.New("metric key is required")
	ErrNilColumn = errors.New("column record is required")
)

// optimizeMask 是位图压缩周期：每新增 4096 个 step 执行一次 RunOptimize，
// 将连续游程收敛为 run 容器。
const optimizeMask = 0xFFF

// Decision 是数据记录的接收结论。
type Decision int8

const (
	// DecisionAccept：新 step，提交阶段将登记位图并更新极值。
	DecisionAccept Decision = iota
	// DecisionSkip：重复或早于等于 min_step 水位线的 step，静默跳过
	//（由调用方输出 debug 日志）。
	DecisionSkip
)

// Precheck 是只读预检的结论，描述一条记录提交后将要发生的状态变更，
// 供调用方在两阶段之间完成容量预留与列记录构建：
//
//	记录产出 = NeedsColumn ? [ColumnRecord, 数据记录?] : [数据记录?]
//
// 注意 NeedsColumn 与 step 结论可以同时成立：即使 step 被拒绝（重复、过期
// 或超限），新列仍会先定义。
type Precheck struct {
	// NeedsColumn: key 尚未定义，提交阶段需要先构建并传入列定义。
	NeedsColumn bool
	// Decision: 数据记录的接纳结论；Err 非空时无意义。
	Decision Decision
	// Err: 确定性错误（wrap ErrStepLimit / ErrTypeMismatch），非 nil 时
	// 该条记录不产生数据记录。
	Err error
}

// RunMetrics 管理一个会话内全部指标的状态。
type RunMetrics struct {
	metrics map[string]*metric
}

// metric 是单个指标 key 的状态。scalar 非空表示标量列，否则为媒体列。
type metric struct {
	column  *columnv1.ColumnRecord
	minStep int64
	steps   *roaring64.Bitmap
	scalar  *scalarState
}

// scalarState 追踪标量极值；仅有限浮点数参与更新（NaN/Inf 忽略）。
type scalarState struct {
	has    bool
	latest float64
	max    float64
	min    float64
}

// New 创建空的指标状态机。
func New() *RunMetrics {
	return &RunMetrics{metrics: make(map[string]*metric)}
}

// Defined 报告 key 是否已定义（供 upsert_columns 的预检使用）。
func (m *RunMetrics) Defined(key string) bool {
	_, ok := m.metrics[key]
	return ok
}

// Len 返回已定义的指标数量。
func (m *RunMetrics) Len() int { return len(m.metrics) }

// MinStep 返回 key 的 step 水位线；key 未定义时第二个返回值为 false。
// 水位线由 Resume 恢复设置，普通自动定义的列为 -1。
func (m *RunMetrics) MinStep(key string) (int64, bool) {
	mt, ok := m.metrics[key]
	if !ok {
		return 0, false
	}
	return mt.minStep, true
}

// DefineScalar 定义一个标量指标。minStep 是 step 下界（含），
// 普通定义为 -1（即允许 step ≥ 0），Resume 恢复时取云端已到达的 step。
func (m *RunMetrics) DefineScalar(key string, col *columnv1.ColumnRecord, minStep int64) error {
	return m.define(key, col, minStep, true)
}

// DefineMedia 定义一个媒体指标。minStep 语义同 DefineScalar。
// 媒体指标不保存存储路径，发送时统一从 RunContext.MediaDir 取。
func (m *RunMetrics) DefineMedia(key string, col *columnv1.ColumnRecord, minStep int64) error {
	return m.define(key, col, minStep, false)
}

func (m *RunMetrics) define(key string, col *columnv1.ColumnRecord, minStep int64, scalar bool) error {
	if key == "" {
		return ErrEmptyKey
	}
	if col == nil {
		return fmt.Errorf("%w: key %q", ErrNilColumn, key)
	}
	if _, exists := m.metrics[key]; exists {
		return fmt.Errorf("%w: %q", ErrColumnExists, key)
	}
	var st *scalarState
	if scalar {
		st = &scalarState{}
	}
	m.metrics[key] = &metric{
		column:  col,
		minStep: minStep,
		steps:   roaring64.New(),
		scalar:  st,
	}
	return nil
}

// ---- 只读预检 ----

// PrecheckScalar 对一条标量记录做只读预检，不产生任何状态变更。
func (m *RunMetrics) PrecheckScalar(rec *datav1.ScalarRecord) Precheck {
	mt, ok := m.metrics[rec.Key]
	if !ok {
		d, err := freshVerdict(rec.Step)
		return Precheck{NeedsColumn: true, Decision: d, Err: err}
	}
	if err := mt.ensureTypeMatch(rec.Type); err != nil {
		return Precheck{Decision: DecisionSkip, Err: err}
	}
	d, err := mt.verdict(rec.Step)
	return Precheck{Decision: d, Err: err}
}

// PrecheckMedia 对一条媒体记录做只读预检，不产生任何状态变更。
func (m *RunMetrics) PrecheckMedia(rec *datav1.MediaRecord) Precheck {
	mt, ok := m.metrics[rec.Key]
	if !ok {
		d, err := freshVerdict(rec.Step)
		return Precheck{NeedsColumn: true, Decision: d, Err: err}
	}
	if err := mt.ensureTypeMatch(rec.Type); err != nil {
		return Precheck{Decision: DecisionSkip, Err: err}
	}
	d, err := mt.verdict(rec.Step)
	return Precheck{Decision: d, Err: err}
}

// freshVerdict 评估自动定义（minStep=-1、空位图）情形下的 step 结论。
func freshVerdict(step int64) (Decision, error) {
	if step <= -1 {
		return DecisionSkip, nil
	}
	if step > MaxStep {
		return DecisionSkip, fmt.Errorf("step %d: %w (max %d)", step, ErrStepLimit, MaxStep)
	}
	return DecisionAccept, nil
}

// ---- 正式提交 ----

// CommitScalar 提交一条标量记录并返回是否接纳了数据记录。
//
// col 仅在该 key 为新列时非 nil（由调用方在预检 NeedsColumn 时用
// builder.BuildAutoColumn 构建）；已定义的 key 传 nil。
// 返回 (false, nil) 表示重复/过期 step 被静默跳过；
// 返回 (false, err) 为确定性错误（类型冲突、step 超限、流程错误）。
func (m *RunMetrics) CommitScalar(rec *datav1.ScalarRecord, col *columnv1.ColumnRecord) (bool, error) {
	mt, err := m.ensureMetric(rec.Key, col, true)
	if err != nil {
		return false, err
	}
	if err = mt.ensureTypeMatch(rec.Type); err != nil {
		return false, err
	}
	accepted, err := mt.accept(rec.Step)
	if err != nil || !accepted {
		return false, err
	}
	// Resume 占位媒体列（UNSPECIFIED）收到标量记录时类型校验放行，
	// 但占位列不更新极值。
	if mt.scalar != nil {
		mt.scalar.update(rec.Key, rec.Value.GetNumber())
	}
	return true, nil
}

// CommitMedia 提交一条媒体记录；col 语义同 CommitScalar。
// 媒体指标没有可更新的状态。
func (m *RunMetrics) CommitMedia(rec *datav1.MediaRecord, col *columnv1.ColumnRecord) (bool, error) {
	mt, err := m.ensureMetric(rec.Key, col, false)
	if err != nil {
		return false, err
	}
	if err = mt.ensureTypeMatch(rec.Type); err != nil {
		return false, err
	}
	return mt.accept(rec.Step)
}

// ensureMetric 取出或定义指标；col 非 nil 且 key 未定义时完成定义。
func (m *RunMetrics) ensureMetric(key string, col *columnv1.ColumnRecord, scalar bool) (*metric, error) {
	if mt, ok := m.metrics[key]; ok {
		if col != nil {
			return nil, fmt.Errorf("%w: %q", ErrColumnExists, key)
		}
		return mt, nil
	}
	if col == nil {
		return nil, fmt.Errorf("%w: %q", ErrColumnMissing, key)
	}
	if err := m.define(key, col, -1, scalar); err != nil {
		return nil, err
	}
	return m.metrics[key], nil
}

// ---- 单指标语义 ----

// ensureTypeMatch 校验列类型与记录类型一致。
//
// 媒体列的例外：Resume 恢复的降级列类型为 UNSPECIFIED（无法得知具体类型），
// 对任意媒体类型放行。
func (mt *metric) ensureTypeMatch(t columnv1.ColumnType) error {
	if mt.scalar == nil && mt.column.ColumnType == columnv1.ColumnType_COLUMN_TYPE_UNSPECIFIED {
		console.Debugf("Media metric '%s' has not been defined, maybe it's defined in resume", mt.column.ColumnKey)
		return nil
	}
	if mt.column.ColumnType != t {
		return fmt.Errorf("%w: metric '%s' has already been defined as %s, not %s",
			ErrTypeMismatch, mt.column.ColumnKey, mt.column.ColumnType, t)
	}
	return nil
}

// verdict 判定顺序：
//  1. step ≤ minStep → Skip（先于超限判定，因此水位线之上的病态 step 也走 Skip）
//  2. step > MaxStep → ErrStepLimit
//  3. 位图已存在 → Skip
//  4. 其余 → Accept
func (mt *metric) verdict(step int64) (Decision, error) {
	if step <= mt.minStep {
		return DecisionSkip, nil
	}
	if step > MaxStep {
		return DecisionSkip, fmt.Errorf("step %d: %w (max %d)", step, ErrStepLimit, MaxStep)
	}
	if mt.steps.Contains(uint64(step)) {
		return DecisionSkip, nil
	}
	return DecisionAccept, nil
}

// accept 登记 step 并按周期压缩位图；false 且 err 为 nil 表示重复/过期，
// false 且 err 非 nil 表示超限。
func (mt *metric) accept(step int64) (bool, error) {
	d, err := mt.verdict(step)
	if err != nil {
		return false, err
	}
	if d != DecisionAccept {
		return false, nil
	}
	mt.steps.Add(uint64(step))
	// 周期性把连续游程压缩为 run 容器，收敛常驻内存
	if mt.steps.GetCardinality()&optimizeMask == 0 {
		mt.steps.RunOptimize()
	}
	return true, nil
}

// update 更新标量极值；仅有限浮点数参与（NaN/±Inf 忽略并输出 debug 日志）。
func (s *scalarState) update(key string, value float64) {
	if !isFinite(value) {
		console.Debugf("Invalid scalar value: %v for metric '%s', ignored when updating.", value, key)
		return
	}
	s.latest = value
	if !s.has || value > s.max {
		s.max = value
	}
	if !s.has || value < s.min {
		s.min = value
	}
	s.has = true
}

func isFinite(v float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0)
}

// SystemKey 判断指标 key 是否为系统指标（前缀 __swanlab__.），
// 转发自 builder 以方便本包的消费方（Resume 摘要解析）。
func SystemKey(key string) bool { return builder.IsSystemKey(key) }
