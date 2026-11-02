package counter

import (
	"sync"
	"testing"
)

func TestCounterBasics(t *testing.T) {
	c := New(0)
	if got := c.Value(); got != 0 {
		t.Fatalf("Value() = %d, want 0", got)
	}
	if got := c.Inc(); got != 1 {
		t.Fatalf("Inc() = %d, want 1", got)
	}
	if got := c.Inc(); got != 2 {
		t.Fatalf("Inc() = %d, want 2", got)
	}
	if got := c.IncBy(10); got != 12 {
		t.Fatalf("IncBy(10) = %d, want 12", got)
	}
	if got := c.Dec(); got != 11 {
		t.Fatalf("Dec() = %d, want 11", got)
	}
	if got := c.DecBy(5); got != 6 {
		t.Fatalf("DecBy(5) = %d, want 6", got)
	}
	c.Reset(100)
	if got := c.Value(); got != 100 {
		t.Fatalf("after Reset(100), Value() = %d, want 100", got)
	}
	if got := c.GetAndReset(-1); got != 100 {
		t.Fatalf("GetAndReset() = %d, want 100", got)
	}
	if got := c.Value(); got != -1 {
		t.Fatalf("after GetAndReset(-1), Value() = %d, want -1", got)
	}
}

func TestCounterZeroValue(t *testing.T) {
	var c Counter
	if got := c.Value(); got != 0 {
		t.Fatalf("零值 Value() = %d, want 0", got)
	}
	if got := c.Inc(); got != 1 {
		t.Fatalf("零值 Inc() = %d, want 1", got)
	}
}

func TestCounterNegativeInitial(t *testing.T) {
	// 离线同步的反向计数从 -100 起递减
	c := New(-100)
	if got := c.Dec(); got != -101 {
		t.Fatalf("Dec() = %d, want -101", got)
	}
}

// TestCounterConcurrentInc 验证并发 Inc 分配出的值两两不同（-race 下检查数据竞争）。
func TestCounterConcurrentInc(t *testing.T) {
	const goroutines, perG = 8, 1000
	c := New(0)

	var wg sync.WaitGroup
	seen := make(chan int64, goroutines*perG)
	for g := 0; g < goroutines; g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := 0; i < perG; i++ {
				seen <- c.Inc()
			}
		}()
	}
	wg.Wait()
	close(seen)

	uniq := make(map[int64]struct{}, goroutines*perG)
	for v := range seen {
		if _, dup := uniq[v]; dup {
			t.Fatalf("并发 Inc 分配了重复值 %d", v)
		}
		uniq[v] = struct{}{}
	}
	if len(uniq) != goroutines*perG {
		t.Fatalf("共分配 %d 个唯一值, want %d", len(uniq), goroutines*perG)
	}
}
