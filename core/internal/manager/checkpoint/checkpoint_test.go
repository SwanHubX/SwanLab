package checkpoint

import (
	"os"
	"path/filepath"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"go.uber.org/goleak"
	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/store"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

func writeTestStore(t *testing.T) (*store.AsyncWriter, string) {
	t.Helper()
	dir := t.TempDir()
	runFile := filepath.Join(dir, "run.swanlab")
	writer, err := store.Create(runFile)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	async, err := store.NewAsyncWriter(writer, store.AsyncOptions{
		Queue:         16,
		FlushBytes:    1,
		FlushInterval: 2 * time.Millisecond,
		IndexStride:   1,
	})
	if err != nil {
		_ = writer.Close()
		t.Fatalf("NewAsyncWriter: %v", err)
	}
	t.Cleanup(func() { _ = async.Close() })
	return async, dir
}

func appendStoreRecords(t *testing.T, async *store.AsyncWriter, recs ...*recordv1.Record) {
	t.Helper()
	for i, rec := range recs {
		payload, err := proto.Marshal(rec)
		if err != nil {
			t.Fatalf("Marshal: %v", err)
		}
		async.In() <- store.Request{
			StoreSeq: int64(i + 1),
			Payload:  payload,
		}
	}
	deadline := time.Now().Add(2 * time.Second)
	for time.Now().Before(deadline) {
		seq, _ := async.Flushed()
		if seq >= int64(len(recs)) {
			return
		}
		time.Sleep(time.Millisecond)
	}
	t.Fatalf("store did not flush %d records", len(recs))
}

func sampleStartRecord() *recordv1.Record {
	return &recordv1.Record{
		Num:        -1,
		RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{}},
	}
}

func sampleScalarRecord(num int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Scalar{
			Scalar: &datav1.ScalarRecord{Key: "loss", Step: num},
		},
	}
}

func TestWriterBoundsByWatermarkAndFlushesFinal(t *testing.T) {
	async, dir := writeTestStore(t)
	runFile := filepath.Join(dir, "run.swanlab")
	cpFile := filepath.Join(dir, "run.checkpoint.json")

	appendStoreRecords(t, async, sampleStartRecord(), sampleScalarRecord(1), sampleScalarRecord(2))

	var watermark atomic.Int64
	watermark.Store(1)
	ticked := make(chan struct{}, 8)

	w := New(Options{
		File:     cpFile,
		RunID:    "run-1",
		RunFile:  runFile,
		Store:    async,
		Interval: 5 * time.Millisecond,
		Watermark: func() int64 {
			return watermark.Load()
		},
		OnTick: func() {
			select {
			case ticked <- struct{}{}:
			default:
			}
		},
	})
	w.Start()
	defer w.Stop()

	select {
	case <-ticked:
	case <-time.After(2 * time.Second):
		t.Fatal("timeout waiting for initial checkpoint tick")
	}

	cp, ok, err := store.LoadCheckpoint(cpFile, "run-1")
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint = (%v, %v), want (nil, true)", err, ok)
	}
	flushedSeq, flushedOff := async.Flushed()
	wantOff, err := store.EndOffset(runFile, async.Index(), 1, flushedOff)
	if err != nil {
		t.Fatalf("EndOffset: %v", err)
	}
	if cp.StoreSeq != 1 || cp.Offset != wantOff {
		t.Fatalf("clamped cp = %+v, want seq=1 off=%d", cp, wantOff)
	}

	watermark.Store(flushedSeq)
	w.Stop()
	w.Flush()

	finalCP, ok, err := store.LoadCheckpoint(cpFile, "run-1")
	if err != nil || !ok {
		t.Fatalf("LoadCheckpoint(final) = (%v, %v), want (nil, true)", err, ok)
	}
	if finalCP.StoreSeq != flushedSeq || finalCP.Offset != flushedOff {
		t.Fatalf("final cp = %+v, want seq=%d off=%d", finalCP, flushedSeq, flushedOff)
	}
}

func TestWriterFailureWarnsAndRecovers(t *testing.T) {
	async, dir := writeTestStore(t)
	cpFile := filepath.Join(dir, "run.checkpoint.json")

	// 用同名非空目录制造 WriteCheckpoint 失败。
	if err := os.MkdirAll(cpFile, 0o755); err != nil {
		t.Fatalf("MkdirAll: %v", err)
	}
	if err := os.WriteFile(filepath.Join(cpFile, "block"), []byte("x"), 0o644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	appendStoreRecords(t, async, sampleStartRecord())
	w := New(Options{
		File:      cpFile,
		RunID:     "run-1",
		RunFile:   filepath.Join(dir, "run.swanlab"),
		Store:     async,
		Watermark: func() int64 { return 1 },
	})

	w.Flush()
	if w.hasWrite || !w.warned {
		t.Fatalf("after failed write: hasWrite=%v warned=%v, want (false, true)", w.hasWrite, w.warned)
	}

	if err := os.RemoveAll(cpFile); err != nil {
		t.Fatalf("RemoveAll: %v", err)
	}
	w.Flush()
	if !w.hasWrite || w.warned {
		t.Fatalf("after recovery: hasWrite=%v warned=%v, want (true, false)", w.hasWrite, w.warned)
	}
	cp, ok, err := store.LoadCheckpoint(cpFile, "run-1")
	if err != nil || !ok || cp.StoreSeq != 1 {
		t.Fatalf("LoadCheckpoint = (%+v, %v, %v), want StoreSeq=1", cp, ok, err)
	}
}

func TestWriterSkipsUnchangedInitialSeqAndNilSafe(t *testing.T) {
	async, dir := writeTestStore(t)
	cpFile := filepath.Join(dir, "run.checkpoint.json")
	appendStoreRecords(t, async, sampleStartRecord())

	ticks := 0
	w := New(Options{
		File:           cpFile,
		RunID:          "run-1",
		RunFile:        filepath.Join(dir, "run.swanlab"),
		Store:          async,
		Watermark:      func() int64 { return 1 },
		InitialSeq:     1,
		InitialWritten: true,
		OnTick:         func() { ticks++ },
	})
	w.Flush()
	if ticks != 0 {
		t.Fatalf("ticks = %d, want 0 when seq matches InitialSeq", ticks)
	}

	// nil & disabled safety
	var nilWriter *Writer
	nilWriter.Start()
	nilWriter.Flush()
	nilWriter.Stop()
	nilWriter.RemoveFile()

	// 缺少 Watermark（非 online 模式无上传水位）不得产生任何 checkpoint。
	noWatermarkCP := filepath.Join(dir, "no-watermark.checkpoint.json")
	disabled := New(Options{
		File:     noWatermarkCP,
		RunID:    "run-1",
		RunFile:  filepath.Join(dir, "run.swanlab"),
		Store:    async,
		Interval: 5 * time.Millisecond,
	})
	disabled.Start()
	disabled.Flush()
	disabled.Stop()
	if disabled.hasWrite {
		t.Fatal("writer without watermark must not record checkpoint writes")
	}
	if _, err := os.Stat(noWatermarkCP); !os.IsNotExist(err) {
		t.Fatalf("writer without watermark must not create checkpoint file, got stat err=%v", err)
	}
}

// TestNewDefaultsToTenSecondInterval 验证默认的定时检查间隔与最小通知间隔。
func TestNewDefaultsToTenSecondInterval(t *testing.T) {
	w := New(Options{})
	if got := w.interval; got != 10*time.Second {
		t.Fatalf("default interval = %s, want 10s", got)
	}
	if got := w.minInterval; got != 2*time.Second {
		t.Fatalf("default minInterval = %s, want 2s", got)
	}
	if got := New(Options{Interval: time.Millisecond}).interval; got != time.Millisecond {
		t.Fatalf("explicit interval = %s, want 1ms", got)
	}
	if got := New(Options{Interval: time.Millisecond}).minInterval; got != time.Millisecond {
		t.Fatalf("clamped minInterval = %s, want 1ms", got)
	}
}

func TestWriterNotifyThrottlesAndFlushSyncsFinal(t *testing.T) {
	async, dir := writeTestStore(t)
	runFile := filepath.Join(dir, "run.swanlab")
	cpFile := filepath.Join(dir, "run.checkpoint.json")
	appendStoreRecords(t, async, sampleStartRecord(), sampleScalarRecord(1), sampleScalarRecord(2), sampleScalarRecord(3))

	var watermark atomic.Int64
	type writeCall struct {
		seq     int64
		durable bool
	}
	var mu sync.Mutex
	var calls []writeCall
	wroteCh := make(chan writeCall, 16)

	w := New(Options{
		File:        cpFile,
		RunID:       "run-1",
		RunFile:     runFile,
		Store:       async,
		Interval:    time.Hour,
		MinInterval: 30 * time.Millisecond,
		Watermark:   watermark.Load,
		writeFn: func(path string, cp store.Checkpoint, durable bool) error {
			var err error
			if durable {
				err = store.WriteCheckpoint(path, cp)
			} else {
				err = store.WriteCheckpointFast(path, cp)
			}
			c := writeCall{seq: cp.StoreSeq, durable: durable}
			mu.Lock()
			calls = append(calls, c)
			mu.Unlock()
			select {
			case wroteCh <- c:
			default:
			}
			return err
		},
	})
	w.Start()
	defer w.Stop()

	// 首次 Notify 立即触发一次运行期非同步写入（durable=false）。
	watermark.Store(2)
	w.Notify()
	select {
	case got := <-wroteCh:
		if got.seq != 2 || got.durable {
			t.Fatalf("first Notify write = %+v, want seq=2 durable=false", got)
		}
	case <-time.After(time.Second):
		t.Fatal("timeout waiting for initial Notify checkpoint write")
	}

	// 窗口内突发多次 Notify 应被合并，并在冷却期结束后写出最新水位。
	watermark.Store(3)
	w.Notify()
	watermark.Store(4)
	w.Notify()
	w.Notify()

	select {
	case got := <-wroteCh:
		if got.seq != 4 || got.durable {
			t.Fatalf("coalesced Notify write = %+v, want seq=4 durable=false", got)
		}
	case <-time.After(time.Second):
		t.Fatal("timeout waiting for coalesced Notify checkpoint write")
	}

	// 停止协程后调用 Flush：即便水位仍是 4，因上次写入是非 fsync 的，Flush 必须补写一次 durable=true。
	w.Stop()
	w.Flush()
	// 水位与持久状态未变时，二次 Flush 应为空操作。
	w.Flush()

	mu.Lock()
	defer mu.Unlock()
	if len(calls) != 3 {
		t.Fatalf("total write calls = %d (%+v), want 3 (2 fast + 1 final sync)", len(calls), calls)
	}
	last := calls[len(calls)-1]
	if last.seq != 4 || !last.durable {
		t.Fatalf("final Flush write = %+v, want seq=4 durable=true", last)
	}
}

func TestWriterAdaptiveBackoffOnSlowWrite(t *testing.T) {
	async, dir := writeTestStore(t)
	runFile := filepath.Join(dir, "run.swanlab")
	cpFile := filepath.Join(dir, "run.checkpoint.json")
	appendStoreRecords(t, async, sampleStartRecord(), sampleScalarRecord(1), sampleScalarRecord(2))

	var watermark atomic.Int64
	var writes atomic.Int32
	wroteAt := make(chan time.Time, 8)

	w := New(Options{
		File:        cpFile,
		RunID:       "run-1",
		RunFile:     runFile,
		Store:       async,
		Interval:    500 * time.Millisecond,
		MinInterval: 10 * time.Millisecond,
		Watermark:   watermark.Load,
		writeFn: func(path string, cp store.Checkpoint, _ bool) error {
			// 模拟慢盘/NFS 耗时 25ms -> 冷却期应自动放大至 >= 100ms (25ms * 4)。
			time.Sleep(25 * time.Millisecond)
			writes.Add(1)
			wroteAt <- time.Now()
			return store.WriteCheckpointFast(path, cp)
		},
	})
	w.Start()
	defer w.Stop()

	watermark.Store(1)
	w.Notify()
	t1 := <-wroteAt

	// 立即推进水位并通知：由于上次写盘耗时 25ms，第二次写盘必须等待 >= ~90ms（受 4x 退避约束）。
	watermark.Store(2)
	w.Notify()
	t2 := <-wroteAt
	gap := t2.Sub(t1)
	if gap < 80*time.Millisecond {
		t.Fatalf("adaptive backoff gap = %v, want >= 80ms after 25ms slow write", gap)
	}
}

func TestWriterStopBeforeStartAndConcurrentLifecycle(t *testing.T) {
	// Stop -> Start -> Stop：先 Stop 后再 Start 与 Notify，不得写盘且不卡住。
	for range 16 {
		async, dir := writeTestStore(t)
		appendStoreRecords(t, async, sampleStartRecord())
		var writes atomic.Int32
		w := New(Options{
			File:      filepath.Join(dir, "run.checkpoint.json"),
			RunID:     "run-1",
			RunFile:   filepath.Join(dir, "run.swanlab"),
			Store:     async,
			Interval:  5 * time.Millisecond,
			Watermark: func() int64 { return 1 },
			writeFn: func(string, store.Checkpoint, bool) error {
				writes.Add(1)
				return nil
			},
		})
		w.Stop()
		w.Start()
		w.Notify()
		time.Sleep(10 * time.Millisecond)
		w.Stop()
		w.Stop()
		if got := writes.Load(); got != 0 {
			t.Fatalf("checkpoint writes = %d after Stop -> Start, want 0", got)
		}
	}

	// 并发 Start / Stop 在超时预算内退出。
	for range 24 {
		async, dir := writeTestStore(t)
		appendStoreRecords(t, async, sampleStartRecord())
		w := New(Options{
			File:      filepath.Join(dir, "run.checkpoint.json"),
			RunID:     "run-1",
			RunFile:   filepath.Join(dir, "run.swanlab"),
			Store:     async,
			Interval:  5 * time.Millisecond,
			Watermark: func() int64 { return 1 },
			writeFn:   func(string, store.Checkpoint, bool) error { return nil },
		})

		var wg sync.WaitGroup
		wg.Add(4)
		go func() { defer wg.Done(); w.Start() }()
		go func() { defer wg.Done(); w.Start() }()
		go func() { defer wg.Done(); w.Stop() }()
		go func() { defer wg.Done(); w.Stop() }()

		done := make(chan struct{})
		go func() {
			wg.Wait()
			close(done)
		}()
		select {
		case <-done:
		case <-time.After(2 * time.Second):
			t.Fatal("concurrent Writer Start/Stop timed out")
		}
		w.Stop()
	}
}
