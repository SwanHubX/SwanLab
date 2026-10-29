package metrics

import (
	"fmt"
	"math/rand"
	"testing"
	"time"

	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
)

// step 去重基准：200 key × 5000 step（100 万 step），场景与规模对齐
// tests/benchmark/sdk/internal/core_python/metrics/bench_metrics_steps.py。
// Steps* 走 CommitScalar 生产路径，Accept* 只测 step 去重（同 try_accept_step）。
// 一次迭代跑完整套 ops，看 ns/step；固定 -benchtime=1x -count=3。
// 运行：go test -bench='Steps|Accept' -benchtime=1x -count=3 ./internal/metrics
const (
	benchKeys  = 200
	benchSteps = 5000
)

func benchDefine(b *testing.B) *RunMetrics {
	m := New()
	for i := 0; i < benchKeys; i++ {
		key := fmt.Sprintf("k%d", i)
		col := &columnv1.ColumnRecord{ColumnKey: key, ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR}
		if err := m.DefineScalar(key, col, -1); err != nil {
			b.Fatal(err)
		}
	}
	return m
}

// genOps 生成与 bench_metrics_steps.py 的 _gen_ops 一致的单个 key 操作序列。
func genOps(workload string, keyIndex int) []int64 {
	switch workload {
	case "sequential":
		out := make([]int64, 0, benchSteps)
		for s := int64(1); s <= benchSteps; s++ {
			out = append(out, s)
		}
		return out
	case "backfill":
		out := make([]int64, 0, benchSteps*6/5)
		for s := int64(1); s <= benchSteps; s++ {
			out = append(out, s)
			if s%10 == 0 {
				out = append(out, s/2) // 重复：必拒绝
				if s > 5 {
					out = append(out, s-5) // 回填：必拒绝
				}
			}
		}
		return out
	default: // sparse
		rng := rand.New(rand.NewSource(int64(keyIndex)))
		seen := make(map[int64]struct{}, benchSteps)
		out := make([]int64, 0, benchSteps)
		for len(out) < benchSteps {
			s := rng.Int63n(MaxStep) + 1
			if _, dup := seen[s]; dup {
				continue
			}
			seen[s] = struct{}{}
			out = append(out, s)
		}
		return out
	}
}

func genAll(workload string) [][]int64 {
	ops := make([][]int64, benchKeys)
	for k := range ops {
		ops[k] = genOps(workload, k)
	}
	return ops
}

// benchRun 走 CommitScalar 完整接收路径。
func benchRun(b *testing.B, m *RunMetrics, ops [][]int64) {
	b.ResetTimer()
	start := time.Now()
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
	elapsed := time.Since(start).Nanoseconds()
	b.StopTimer()
	total := accepted + rejected
	if total > 0 {
		// b.N 次迭代各跑完整套 ops，ns/op 因此是整套耗时；这里按 step 折算。
		b.ReportMetric(float64(elapsed)/float64(total), "ns/step")
	}
	b.ReportMetric(float64(total), "ops")
	b.ReportMetric(float64(rejected), "rejected")
}

// benchRunAccept 只计时 step 去重，对齐 Python 的 try_accept_step。
func benchRunAccept(b *testing.B, m *RunMetrics, ops [][]int64) {
	b.ResetTimer()
	start := time.Now()
	accepted, total := 0, 0
	for k, steps := range ops {
		mt := m.metrics[fmt.Sprintf("k%d", k)]
		for _, s := range steps {
			total++
			if ok, err := mt.accept(s); err == nil && ok {
				accepted++
			}
		}
	}
	elapsed := time.Since(start).Nanoseconds()
	b.StopTimer()
	if total > 0 {
		b.ReportMetric(float64(elapsed)/float64(total), "ns/step")
	}
	// bytes/key：位图原生序列化体积（与 pyroaring portable 格式不同，仅同实现内对比）
	var size uint64
	for _, mt := range m.metrics {
		size += mt.steps.GetSerializedSizeInBytes()
	}
	b.ReportMetric(float64(size)/float64(benchKeys), "bytes/key")
	b.ReportMetric(float64(total), "ops")
	b.ReportMetric(float64(total-accepted), "rejected")
}

// BenchmarkStepsSequential: 每个 key 顺序写入 1..N（真实训练形态）。
func BenchmarkStepsSequential(b *testing.B) { benchRun(b, benchDefine(b), genAll("sequential")) }

// BenchmarkStepsBackfill: 顺序写入 + 周期性重复/回填探测（去重拒绝路径）。
func BenchmarkStepsBackfill(b *testing.B) { benchRun(b, benchDefine(b), genAll("backfill")) }

// BenchmarkStepsSparse: 每个 key 在 [1, 2^48) 内随机不重复 step（病态稀疏）。
func BenchmarkStepsSparse(b *testing.B) { benchRun(b, benchDefine(b), genAll("sparse")) }

// BenchmarkAcceptSequential: 同 StepsSequential，只测 step 去重（对齐 try_accept_step）。
func BenchmarkAcceptSequential(b *testing.B) { benchRunAccept(b, benchDefine(b), genAll("sequential")) }

// BenchmarkAcceptBackfill: 同 StepsBackfill，只测 step 去重。
func BenchmarkAcceptBackfill(b *testing.B) { benchRunAccept(b, benchDefine(b), genAll("backfill")) }

// BenchmarkAcceptSparse: 同 StepsSparse，只测 step 去重。
func BenchmarkAcceptSparse(b *testing.B) { benchRunAccept(b, benchDefine(b), genAll("sparse")) }
