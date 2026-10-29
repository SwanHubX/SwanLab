package uploader

import (
	"math"
	"sync"
	"testing"
	"time"

	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

// fakeClock 是测试用可推进时钟，避免真实等待。
type fakeClock struct{ t time.Time }

func newFakeClock() *fakeClock {
	return &fakeClock{t: time.Unix(1_700_000_000, 0)}
}

func (c *fakeClock) now() time.Time { return c.t }

func (c *fakeClock) advance(d time.Duration) { c.t = c.t.Add(d) }

func almostEqual(a, b float64) bool { return math.Abs(a-b) < 1e-6 }

func TestTrackerRecords(t *testing.T) {
	tr := NewTracker()
	tr.AddTotalRecords(3)
	tr.AddTotalRecords(0)
	tr.AddTotalRecords(-2)
	tr.AdvanceRecords(2)
	tr.AdvanceRecords(-1)

	stats := tr.Snapshot()
	if stats.TotalRecords != 3 || stats.UploadedRecords != 2 {
		t.Fatalf("records = %d/%d, want 3/2", stats.UploadedRecords, stats.TotalRecords)
	}
	if stats.State != operationv1.CoreState_CORE_STATE_NOT_STARTED {
		t.Fatalf("initial state = %v", stats.State)
	}

	tr.SetState(operationv1.CoreState_CORE_STATE_FINISHED)
	if got := tr.Snapshot().State; got != operationv1.CoreState_CORE_STATE_FINISHED {
		t.Fatalf("state = %v, want FINISHED", got)
	}
}

func TestTrackerFileLifecycle(t *testing.T) {
	clock := newFakeClock()
	tr := NewTracker()
	tr.now = clock.now

	tr.AddFile("f1", "/tmp/a", 1000)
	tr.AddFile("f1", "/tmp/a", 1000) // 幂等
	tr.AddFile("f2", "/tmp/zero", 0) // total<=0 忽略

	stats := tr.Snapshot()
	if stats.TotalNumber != 1 || stats.TotalSize != 1000 || stats.UploadedNumber != 0 {
		t.Fatalf("after add: number=%d size=%d uploaded=%d", stats.TotalNumber, stats.TotalSize, stats.UploadedNumber)
	}

	tr.UpdateFileProgress("f1", "part-1", 300)
	tr.UpdateFileProgress("f1", "part-1", 300) // 重复上报同一分片不累计
	tr.UpdateFileProgress("f1", "part-2", 200)
	tr.UpdateFileProgress("f1", "part-1", 100)      // 负增量且 current>0，忽略
	tr.UpdateFileProgress("missing", "part-1", 100) // 未知文件，忽略

	if got := tr.Snapshot().UploadedSize; got != 500 {
		t.Fatalf("uploaded size = %d, want 500", got)
	}

	tr.FinishFile("f1")
	stats = tr.Snapshot()
	if stats.UploadedNumber != 1 || stats.UploadedSize != 1000 {
		t.Fatalf("after finish: number=%d size=%d, want 1/1000", stats.UploadedNumber, stats.UploadedSize)
	}
	if len(stats.Files) != 0 {
		t.Fatalf("finished file should not be visible: %d", len(stats.Files))
	}

	tr.FinishFile("f1") // 幂等
	if got := tr.Snapshot().UploadedNumber; got != 1 {
		t.Fatalf("uploaded number = %d, want 1", got)
	}
}

func TestTrackerFailFile(t *testing.T) {
	tr := NewTracker()
	tr.AddFile("f1", "/tmp/a", 100)
	tr.UpdateFileProgress("f1", "part-1", 40)
	tr.FailFile("f1")

	stats := tr.Snapshot()
	if stats.UploadedNumber != 0 || stats.UploadedSize != 40 {
		t.Fatalf("failed file: number=%d size=%d", stats.UploadedNumber, stats.UploadedSize)
	}
	if len(stats.Files) != 0 {
		t.Fatalf("failed file should be removed: %d", len(stats.Files))
	}
}

func TestTrackerVisibleFiles(t *testing.T) {
	tr := NewTracker()
	for _, key := range []string{"f1", "f2", "f3", "f4", "f5", "f6", "f7"} {
		tr.AddFile(key, "/tmp/"+key, 100)
	}
	tr.UpdateFileProgress("f2", "part", 10)
	tr.FinishFile("f1")

	stats := tr.Snapshot()
	if len(stats.Files) != defaultVisibleFiles {
		t.Fatalf("visible files = %d, want %d", len(stats.Files), defaultVisibleFiles)
	}
	wantPaths := []string{"/tmp/f2", "/tmp/f3", "/tmp/f4", "/tmp/f5", "/tmp/f6"}
	for i, path := range wantPaths {
		if stats.Files[i].Path != path {
			t.Errorf("visible[%d] = %s, want %s", i, stats.Files[i].Path, path)
		}
	}

	tr.limit = 0
	if got := len(tr.Snapshot().Files); got != 0 {
		t.Fatalf("limit=0 should hide files, got %d", got)
	}
}

func TestTrackerRateEMA(t *testing.T) {
	clock := newFakeClock()
	tr := NewTracker()
	tr.now = clock.now
	tr.AddFile("f1", "/tmp/a", 10_000)

	// 第一个窗口：0.6s 内 600 字节 → 1000 B/s。
	clock.advance(600 * time.Millisecond)
	tr.UpdateFileProgress("f1", "part", 600)
	if got := tr.Snapshot().Rate; !almostEqual(got, 1000) {
		t.Fatalf("first window rate = %v, want 1000", got)
	}

	// 窗口未满：速率不变。
	clock.advance(400 * time.Millisecond)
	tr.UpdateFileProgress("f1", "part", 700)
	if got := tr.Snapshot().Rate; !almostEqual(got, 1000) {
		t.Fatalf("rate before window boundary = %v, want 1000", got)
	}

	// 第二个窗口：0.6s 内 400 字节 → 瞬时 666.67，EMA = 0.3*666.67 + 0.7*1000 = 900。
	clock.advance(200 * time.Millisecond)
	tr.UpdateFileProgress("f1", "part", 1000)
	stats := tr.Snapshot()
	if got := stats.Rate; !almostEqual(got, 900) {
		t.Fatalf("EMA rate = %v, want 900", got)
	}
	if len(stats.Files) != 1 || !almostEqual(stats.Files[0].Rate, 900) {
		t.Fatalf("file rate = %+v, want 900", stats.Files)
	}

	// 一个窗口无进展后按 (1-alpha) 衰减。
	clock.advance(trackerRateWindow)
	if got := tr.Snapshot().Rate; !almostEqual(got, 630) {
		t.Fatalf("decayed rate = %v, want 630", got)
	}
	clock.advance(trackerRateWindow)
	if got := tr.Snapshot().Rate; !almostEqual(got, 441) {
		t.Fatalf("decayed rate = %v, want 441", got)
	}
}

func TestTrackerConcurrent(t *testing.T) {
	tr := NewTracker()
	var wg sync.WaitGroup
	for i := 0; i < 8; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := 0; j < 1000; j++ {
				tr.AddTotalRecords(1)
				tr.AdvanceRecords(1)
			}
		}()
	}
	wg.Wait()

	stats := tr.Snapshot()
	if stats.TotalRecords != 8000 || stats.UploadedRecords != 8000 {
		t.Fatalf("records = %d/%d, want 8000/8000", stats.UploadedRecords, stats.TotalRecords)
	}
}
