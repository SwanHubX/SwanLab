package metrics

import (
	"fmt"
	"math/rand"
	"testing"

	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
)

// step 去重基准：三种工作负载 × 200 key，每 key 500 步，单负载 10 万 step。
// 走 CommitScalar 完整路径（类型校验 + 位图 + 极值），打印 ns/op。
// 运行：go test -bench=BenchSteps -benchtime=1x ./internal/metrics
const (
	benchKeys  = 200
	benchSteps = 500
)

func benchDefine(b *testing.B) *RunMetrics {
	m := New()
	for i := 0; i < benchKeys; i++ {
		if err := m.DefineScalar(fmt.Sprintf("k%d", i),
			&columnv1.ColumnRecord{ColumnKey: fmt.Sprintf("k%d", i), ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR}, -1); err != nil {
			b.Fatal(err)
		}
	}
	return m
}

func benchRun(b *testing.B, m *RunMetrics, ops [][]int64) {
	b.ResetTimer()
	accepted, rejected := 0, 0
	for k, steps := range ops {
		key := fmt.Sprintf("k%d", k)
		for _, s := range steps {
			if ok, err := m.CommitScalar(&datav1.ScalarRecord{
				Key: key, Step: s, Type: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
				Value: &datav1.ScalarValue{Number: float64(s)},
			}, nil); err == nil && ok {
				accepted++
			} else {
				rejected++
			}
		}
	}
	b.StopTimer()
	b.ReportMetric(float64(accepted+rejected), "ops")
	b.ReportMetric(float64(rejected), "rejected")
}

// BenchmarkStepsSequential: 每个 key 顺序写入 1..N（真实训练形态）。
func BenchmarkStepsSequential(b *testing.B) {
	m := benchDefine(b)
	ops := make([][]int64, benchKeys)
	for k := range ops {
		ops[k] = make([]int64, benchSteps)
		for s := int64(0); s < benchSteps; s++ {
			ops[k][s] = s + 1 // 顺序 1..N
		}
	}
	benchRun(b, m, ops)
}

// BenchmarkStepsBackfill: 顺序写入 + 周期性重复/回填探测（去重拒绝路径）。
func BenchmarkStepsBackfill(b *testing.B) {
	m := benchDefine(b)
	ops := make([][]int64, benchKeys)
	for k := range ops {
		for s := int64(1); s <= benchSteps; s++ {
			ops[k] = append(ops[k], s)
			if s%10 == 0 {
				ops[k] = append(ops[k], s/2) // 重复：必拒绝
				if s > 5 {
					ops[k] = append(ops[k], s-5) // 回填：必拒绝
				}
			}
		}
	}
	benchRun(b, m, ops)
}

// BenchmarkStepsSparse: 每个 key 在 [1, 2^48) 内随机不重复 step（病态稀疏）。
func BenchmarkStepsSparse(b *testing.B) {
	m := benchDefine(b)
	ops := make([][]int64, benchKeys)
	for k := range ops {
		rng := rand.New(rand.NewSource(int64(k)))
		seen := make(map[int64]struct{}, benchSteps)
		steps := make([]int64, 0, benchSteps)
		for len(steps) < benchSteps {
			s := rng.Int63n(MaxStep) + 1
			if _, dup := seen[s]; dup {
				continue
			}
			seen[s] = struct{}{}
			steps = append(steps, s)
		}
		ops[k] = steps
	}
	benchRun(b, m, ops)
}
