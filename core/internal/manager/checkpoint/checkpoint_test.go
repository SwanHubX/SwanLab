package checkpoint

import (
	"os"
	"path/filepath"
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
			StoreSeq: uint64(i + 1),
			Payload:  payload,
		}
	}
	deadline := time.Now().Add(2 * time.Second)
	for time.Now().Before(deadline) {
		seq, _ := async.Flushed()
		if seq >= uint64(len(recs)) {
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

	var watermark atomic.Uint64
	watermark.Store(1)
	ticked := make(chan struct{}, 8)

	w := New(Options{
		File:     cpFile,
		RunID:    "run-1",
		RunFile:  runFile,
		Store:    async,
		Interval: 5 * time.Millisecond,
		Watermark: func() uint64 {
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
		Watermark: func() uint64 { return 1 },
	})

	w.Persist()
	if w.hasWrite || !w.warned {
		t.Fatalf("after failed write: hasWrite=%v warned=%v, want (false, true)", w.hasWrite, w.warned)
	}

	if err := os.RemoveAll(cpFile); err != nil {
		t.Fatalf("RemoveAll: %v", err)
	}
	w.Persist()
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
		Watermark:      func() uint64 { return 1 },
		InitialSeq:     1,
		InitialWritten: true,
		OnTick:         func() { ticks++ },
	})
	w.Persist()
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
	disabled.Persist()
	disabled.Flush()
	disabled.Stop()
	if disabled.hasWrite {
		t.Fatal("writer without watermark must not record checkpoint writes")
	}
	if _, err := os.Stat(noWatermarkCP); !os.IsNotExist(err) {
		t.Fatalf("writer without watermark must not create checkpoint file, got stat err=%v", err)
	}
}
