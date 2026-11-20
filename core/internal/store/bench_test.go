package store

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"testing"
	"time"
)

// 本地写盘基准：100 万条 120 字节记录，规模与 payload 大小对齐
// tests/benchmark/python_core/store/bench_store_fsync.py。
//
//   - BenchmarkStoreWrite 走同步路径（Create + WriteRecord + Flush + Close）；
//   - BenchmarkStoreRead  回读 BenchmarkStoreWrite 生成的同构文件；
//   - BenchmarkStoreAsync 走生产路径（AsyncWriter 有界队列 + 定量刷盘），
//     队列容量分别为默认的 8×record_batch=80,000 与压测档 8,192。
//
// 固定 -benchtime=1x -count=3，取三次中位数。
// 运行：go test -run '^$' -bench 'Store' -benchtime=1x -count=3 ./internal/store
const (
	benchRecords     = 1_000_000
	benchPayloadSize = 120
)

// benchPayload 生成一条可复用的定长 payload（内容无关，长度决定分帧）。
func benchPayload() []byte {
	p := make([]byte, benchPayloadSize)
	for i := range p {
		p[i] = byte(i * 7)
	}
	return p
}

// benchPercentileUs 返回已排序耗时切片的百分位值（微秒）。
func benchPercentileUs(sorted []time.Duration, pct float64) float64 {
	if len(sorted) == 0 {
		return 0
	}
	idx := int(float64(len(sorted)) * pct / 100)
	if idx >= len(sorted) {
		idx = len(sorted) - 1
	}
	return float64(sorted[idx].Nanoseconds()) / 1e3
}

// benchWriteFile 同步写入 benchRecords 条记录并逐条计时，返回每条耗时。
func benchWriteFile(path string, payload []byte) ([]time.Duration, time.Duration, int64, error) {
	w, err := Create(path)
	if err != nil {
		return nil, 0, 0, err
	}
	latencies := make([]time.Duration, 0, benchRecords)
	start := time.Now()
	for n := 0; n < benchRecords; n++ {
		t0 := time.Now()
		if _, err = w.WriteRecord(payload); err != nil {
			w.Close()
			return nil, 0, 0, err
		}
		latencies = append(latencies, time.Since(t0))
	}
	elapsed := time.Since(start)
	if err = w.Flush(); err != nil {
		w.Close()
		return nil, 0, 0, err
	}
	if err = w.Close(); err != nil {
		return nil, 0, 0, err
	}
	info, err := os.Stat(path)
	if err != nil {
		return nil, 0, 0, err
	}
	return latencies, elapsed, info.Size(), nil
}

// BenchmarkStoreWrite 同步写入 100 万条 120 字节记录。
//
// 注意：单条耗时为逐条 time.Now() 采样，含计时开销（本机约 25ns/次），
// 因此 p50 会偏高；吞吐与 ns/rec 以整段 elapsed 为准。
func BenchmarkStoreWrite(b *testing.B) {
	payload := benchPayload()
	dir := b.TempDir()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		path := filepath.Join(dir, fmt.Sprintf("write-%d.swanlab", i))
		latencies, elapsed, size, err := benchWriteFile(path, payload)
		if err != nil {
			b.Fatal(err)
		}
		sort.Slice(latencies, func(a, c int) bool { return latencies[a] < latencies[c] })
		b.StopTimer()
		b.ReportMetric(float64(benchRecords)/elapsed.Seconds(), "rec/s")
		b.ReportMetric(float64(elapsed.Nanoseconds())/benchRecords, "ns/rec")
		b.ReportMetric(float64(size)/benchRecords, "bytes/rec")
		b.ReportMetric(benchPercentileUs(latencies, 50), "p50_us")
		b.ReportMetric(benchPercentileUs(latencies, 95), "p95_us")
		b.ReportMetric(benchPercentileUs(latencies, 99), "p99_us")
		b.StartTimer()
	}
}

// BenchmarkStoreRead 顺序回读同一份 100 万条文件。
//
// 数据准备（写入）在计时之外，只测 Reader.Scan 的解析成本；
// 回读条数与字节数在每次迭代内校验，作为数据完整性断言。
func BenchmarkStoreRead(b *testing.B) {
	payload := benchPayload()
	path := filepath.Join(b.TempDir(), "read.swanlab")
	b.StopTimer()
	if _, _, size, err := benchWriteFile(path, payload); err != nil {
		b.Fatal(err)
	} else {
		b.Logf("prepared %d records, %d bytes", benchRecords, size)
	}
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		r, err := Open(path)
		if err != nil {
			b.Fatal(err)
		}
		count, bytes := 0, int64(0)
		start := time.Now()
		for {
			record, err := r.Scan()
			if err == io.EOF {
				break
			}
			if err != nil {
				r.Close()
				b.Fatalf("scan after %d records: %v", count, err)
			}
			count++
			bytes += int64(len(record))
		}
		elapsed := time.Since(start)
		if err := r.Close(); err != nil {
			b.Fatal(err)
		}
		b.StopTimer()
		if count != benchRecords {
			b.Fatalf("read %d records, want %d", count, benchRecords)
		}
		b.ReportMetric(float64(count)/elapsed.Seconds(), "rec/s")
		b.ReportMetric(float64(elapsed.Nanoseconds())/float64(count), "ns/rec")
		b.ReportMetric(float64(bytes)/float64(count), "bytes/rec")
		b.ReportMetric(float64(bytes)/elapsed.Seconds()/(1<<20), "MiB/s")
		b.StartTimer()
	}
}

// BenchmarkStoreAsync 走生产路径：有界队列 + 后台写协程 + 定量刷盘。
//
// 生产方按 intake 的方式「先试非阻塞发送，队列满再阻塞等待」，只在阻塞时计时，
// 因此 p99/max 反映的是队列满时的入队等待，而不是平常的 channel 发送成本。
// 每条记录单独分配 payload（生产也是逐条新建），堆增量因此包含记录本体的分配。
func BenchmarkStoreAsync(b *testing.B) {
	cases := []struct {
		name  string
		queue int
	}{
		{"queue-80000", 80_000}, // 8 × record_batch(10000)，默认档
		{"queue-8192", 8_192},   // 压测档：制造队列打满
	}
	template := benchPayload()
	for _, tc := range cases {
		b.Run(tc.name, func(b *testing.B) {
			dir := b.TempDir()
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				path := filepath.Join(dir, fmt.Sprintf("async-%d.swanlab", i))
				w, err := Create(path)
				if err != nil {
					b.Fatal(err)
				}
				aw, err := NewAsyncWriter(w, AsyncOptions{Queue: tc.queue})
				if err != nil {
					b.Fatal(err)
				}

				var before, during, after runtime.MemStats
				runtime.GC()
				runtime.ReadMemStats(&before)

				waits := make([]time.Duration, 0, 1024)
				blocked := 0
				start := time.Now()
				for n := 0; n < benchRecords; n++ {
					payload := make([]byte, benchPayloadSize)
					copy(payload, template)
					req := Request{StoreSeq: uint64(n + 1), Payload: payload}
					select {
					case aw.In() <- req:
					default:
						t0 := time.Now()
						aw.In() <- req
						waits = append(waits, time.Since(t0))
						blocked++
					}
				}
				enqueue := time.Since(start)
				// 入队刚结束时队列最满，此刻的 HeapAlloc 近似「队列占用的堆」。
				runtime.ReadMemStats(&during)
				if closeErr := aw.Close(); closeErr != nil {
					b.Fatal(closeErr)
				}
				elapsed := time.Since(start)

				runtime.GC()
				runtime.ReadMemStats(&after)

				info, err := os.Stat(path)
				if err != nil {
					b.Fatal(err)
				}
				sort.Slice(waits, func(a, c int) bool { return waits[a] < waits[c] })
				b.StopTimer()
				b.ReportMetric(float64(benchRecords)/elapsed.Seconds(), "rec/s")
				b.ReportMetric(float64(elapsed.Nanoseconds())/benchRecords, "ns/rec")
				b.ReportMetric(float64(enqueue.Nanoseconds())/benchRecords, "enqueue_ns/rec")
				b.ReportMetric(float64(blocked), "blocked_sends")
				b.ReportMetric(benchPercentileUs(waits, 50), "wait_p50_us")
				b.ReportMetric(benchPercentileUs(waits, 99), "wait_p99_us")
				b.ReportMetric(benchPercentileUs(waits, 100), "wait_max_us")
				b.ReportMetric(float64(info.Size())/benchRecords, "bytes/rec")
				b.ReportMetric(float64(int64(during.HeapAlloc)-int64(before.HeapAlloc))/(1<<20), "heap_full_MiB")
				b.ReportMetric(float64(int64(after.HeapAlloc)-int64(before.HeapAlloc))/(1<<20), "heap_final_MiB")
				b.ReportMetric(float64(after.TotalAlloc-before.TotalAlloc)/benchRecords, "alloc_B/rec")
				b.StartTimer()
			}
		})
	}
}
