package budget

import (
	"context"
	"errors"
	"sync"
	"testing"
	"time"
)

func TestNewNormalizesLimits(t *testing.T) {
	b := New(Limits{Records: -1, Bytes: -5})
	if got := b.Limits(); got != (Limits{}) {
		t.Fatalf("Limits() = %+v, want zero value", got)
	}
	if !b.Limits().Unlimited() {
		t.Fatal("negative limits should be treated as unlimited")
	}
	if _, bytes := b.Usage(); bytes != 0 {
		t.Fatalf("initial bytes = %d, want 0", bytes)
	}
}

func TestTryAcquire(t *testing.T) {
	tests := []struct {
		name    string
		limits  Limits
		first   [2]int64 // records, bytes
		second  [2]int64
		wantOK  bool
		wantRec int64
		wantByt int64
	}{
		{
			name:    "unlimited always fits",
			limits:  Limits{},
			first:   [2]int64{100, 1 << 20},
			second:  [2]int64{100, 1 << 20},
			wantOK:  true,
			wantRec: 200,
			wantByt: 2 << 20,
		},
		{
			name:    "records limit",
			limits:  Limits{Records: 3},
			first:   [2]int64{2, 0},
			second:  [2]int64{2, 0},
			wantOK:  false,
			wantRec: 2,
		},
		{
			name:    "bytes limit",
			limits:  Limits{Bytes: 10},
			first:   [2]int64{1, 6},
			second:  [2]int64{1, 5},
			wantOK:  false,
			wantRec: 1,
			wantByt: 6,
		},
		{
			name:    "exact fit",
			limits:  Limits{Records: 2, Bytes: 10},
			first:   [2]int64{2, 10},
			second:  [2]int64{1, 0},
			wantOK:  false,
			wantRec: 2,
			wantByt: 10,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			b := New(tc.limits)
			if !b.TryAcquire(tc.first[0], tc.first[1]) {
				t.Fatalf("first TryAcquire(%d, %d) = false, want true", tc.first[0], tc.first[1])
			}
			if got := b.TryAcquire(tc.second[0], tc.second[1]); got != tc.wantOK {
				t.Fatalf("second TryAcquire(%d, %d) = %v, want %v", tc.second[0], tc.second[1], got, tc.wantOK)
			}
			records, bytes := b.Usage()
			if records != tc.wantRec || bytes != tc.wantByt {
				t.Fatalf("Usage() = (%d, %d), want (%d, %d)", records, bytes, tc.wantRec, tc.wantByt)
			}
		})
	}
}

func TestTryAcquireTooLarge(t *testing.T) {
	b := New(Limits{Records: 1})
	if b.TryAcquire(2, 0) {
		t.Fatal("TryAcquire(2, 0) = true with limit 1")
	}
	if records, bytes := b.Usage(); records != 0 || bytes != 0 {
		t.Fatalf("Usage() = (%d, %d), want (0, 0)", records, bytes)
	}
}

func TestAcquireWaitsForRelease(t *testing.T) {
	b := New(Limits{Records: 1})
	if !b.TryAcquire(1, 0) {
		t.Fatal("TryAcquire(1, 0) = false, want true")
	}

	acquired := make(chan error, 1)
	go func() { acquired <- b.Acquire(context.Background(), 1, 0) }()

	select {
	case err := <-acquired:
		t.Fatalf("Acquire returned before Release: %v", err)
	case <-time.After(20 * time.Millisecond):
	}

	b.Release(1, 0)
	select {
	case err := <-acquired:
		if err != nil {
			t.Fatalf("Acquire after Release = %v, want nil", err)
		}
	case <-time.After(time.Second):
		t.Fatal("Acquire did not wake after Release")
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("Usage records = %d, want 1", records)
	}
}

func TestAcquireContextCanceled(t *testing.T) {
	b := New(Limits{Records: 1})
	if !b.TryAcquire(1, 0) {
		t.Fatal("TryAcquire(1, 0) = false, want true")
	}
	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
	defer cancel()

	err := b.Acquire(ctx, 1, 0)
	if !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("Acquire err = %v, want context.DeadlineExceeded", err)
	}
	if records, _ := b.Usage(); records != 1 {
		t.Fatalf("Usage records = %d, want 1 (failed acquire must not change state)", records)
	}
}

func TestAcquireTooLarge(t *testing.T) {
	b := New(Limits{Records: 2})
	err := b.Acquire(context.Background(), 3, 0)
	if !errors.Is(err, ErrTooLarge) {
		t.Fatalf("Acquire err = %v, want ErrTooLarge", err)
	}
	if records, _ := b.Usage(); records != 0 {
		t.Fatalf("Usage records = %d, want 0", records)
	}
}

func TestReleasePanicsOnOverRelease(t *testing.T) {
	tests := []struct {
		name      string
		acquire   [2]int64
		release   [2]int64
		wantPanic bool
	}{
		{name: "exact", acquire: [2]int64{2, 10}, release: [2]int64{2, 10}},
		{name: "over records", acquire: [2]int64{2, 10}, release: [2]int64{3, 10}, wantPanic: true},
		{name: "over bytes", acquire: [2]int64{2, 10}, release: [2]int64{2, 11}, wantPanic: true},
		{name: "negative", acquire: [2]int64{2, 10}, release: [2]int64{-1, 0}, wantPanic: true},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			b := New(Limits{})
			b.TryAcquire(tc.acquire[0], tc.acquire[1])
			defer func() {
				got := recover() != nil
				if got != tc.wantPanic {
					t.Fatalf("panic = %v, want %v", got, tc.wantPanic)
				}
			}()
			b.Release(tc.release[0], tc.release[1])
		})
	}
}

func TestAdjust(t *testing.T) {
	b := New(Limits{Records: 10, Bytes: 100})
	if !b.TryAcquire(2, 50) {
		t.Fatal("TryAcquire(2, 50) = false, want true")
	}
	// 实际比估算小：退还差额。
	b.Adjust(1, -20)
	if records, bytes := b.Usage(); records != 3 || bytes != 30 {
		t.Fatalf("Usage after shrink = (%d, %d), want (3, 30)", records, bytes)
	}
	// 实际比估算大：补扣差额。
	b.Adjust(1, 20)
	if records, bytes := b.Usage(); records != 4 || bytes != 50 {
		t.Fatalf("Usage after grow = (%d, %d), want (4, 50)", records, bytes)
	}
}

func TestAdjustPanicsOnOverRefund(t *testing.T) {
	b := New(Limits{})
	b.TryAcquire(1, 10)
	defer func() {
		if recover() == nil {
			t.Fatal("Adjust with excessive refund did not panic")
		}
	}()
	b.Adjust(0, -11)
}

func TestWatermark(t *testing.T) {
	b := New(Limits{Records: 4, Bytes: 100})
	if got := b.Ratio(); got != 0 {
		t.Fatalf("empty Ratio = %v, want 0", got)
	}
	b.TryAcquire(2, 50)
	if got := b.Ratio(); got != 0.5 {
		t.Fatalf("Ratio = %v, want 0.5", got)
	}
	if b.Full() {
		t.Fatal("Full() = true at 50%, want false")
	}
	b.TryAcquire(2, 50)
	if got := b.Ratio(); got != 1 {
		t.Fatalf("Ratio = %v, want 1", got)
	}
	if !b.Full() {
		t.Fatal("Full() = false at 100%, want true")
	}

	unlimited := New(Limits{})
	unlimited.TryAcquire(1<<20, 1<<30)
	if got := unlimited.Ratio(); got != 0 {
		t.Fatalf("unlimited Ratio = %v, want 0", got)
	}
	if unlimited.Full() {
		t.Fatal("unlimited Full() = true, want false")
	}
}

// TestQuotaConservation 验证配额在「预留 → 组件间流转 → 覆盖/终结」全过程中守恒：
// 流转不释放，只有覆盖与终结才归还。
func TestQuotaConservation(t *testing.T) {
	tests := []struct {
		name    string
		actions func(t *testing.T, b *Budget)
		wantRec int64
		wantByt int64
	}{
		{
			name: "transfer keeps quota",
			actions: func(t *testing.T, b *Budget) {
				// intake → upload channel → batcher: 只是换持有者。
				if !b.TryAcquire(2, 20) {
					t.Fatal("TryAcquire = false")
				}
			},
			wantRec: 2,
			wantByt: 20,
		},
		{
			name: "overwrite releases only the replaced record",
			actions: func(t *testing.T, b *Budget) {
				if !b.TryAcquire(2, 20) {
					t.Fatal("TryAcquire = false")
				}
				// 同号覆盖：新记录占 1 条 10 字节，被覆盖的旧记录归还 1 条 10 字节。
				b.Release(1, 10)
				if !b.TryAcquire(1, 10) {
					t.Fatal("TryAcquire for replacement = false")
				}
			},
			wantRec: 2,
			wantByt: 20,
		},
		{
			name: "confirmed releases",
			actions: func(t *testing.T, b *Budget) {
				if !b.TryAcquire(2, 20) {
					t.Fatal("TryAcquire = false")
				}
				b.Release(1, 10)
				b.Release(1, 10)
			},
			wantRec: 0,
			wantByt: 0,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			b := New(Limits{Records: 10, Bytes: 100})
			tc.actions(t, b)
			records, bytes := b.Usage()
			if records != tc.wantRec || bytes != tc.wantByt {
				t.Fatalf("Usage() = (%d, %d), want (%d, %d)", records, bytes, tc.wantRec, tc.wantByt)
			}
		})
	}
}

func TestConcurrentAcquireRelease(t *testing.T) {
	const (
		workers = 8
		rounds  = 200
	)
	b := New(Limits{Records: 4, Bytes: 40})
	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := 0; j < rounds; j++ {
				if err := b.Acquire(context.Background(), 1, 10); err != nil {
					t.Errorf("Acquire = %v", err)
					return
				}
				b.Release(1, 10)
			}
		}()
	}
	wg.Wait()
	if records, bytes := b.Usage(); records != 0 || bytes != 0 {
		t.Fatalf("Usage after all workers = (%d, %d), want (0, 0)", records, bytes)
	}
}
