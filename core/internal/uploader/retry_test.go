package uploader

import (
	"testing"
	"time"
)

// fixedRand 返回恒定随机源，让抖动系数可以手算。
func fixedRand(v float64) func() float64 {
	return func() float64 { return v }
}

func TestDefaultRetryPolicy(t *testing.T) {
	p := DefaultRetryPolicy()
	if p.Base != time.Second || p.Max != 10*time.Second || p.Multiplier != 2 {
		t.Fatalf("unexpected backoff shape: %+v", p)
	}
	if p.JitterLow != 0.8 || p.JitterHigh != 1.2 {
		t.Fatalf("unexpected jitter: %v..%v", p.JitterLow, p.JitterHigh)
	}
	if p.Rand == nil {
		t.Fatal("default policy must provide a random source")
	}
}

func TestRetryPolicyDelaySequence(t *testing.T) {
	p := DefaultRetryPolicy()
	p.JitterLow, p.JitterHigh = 1, 1 // 关闭抖动，只看指数项
	p.Rand = fixedRand(0)

	want := []time.Duration{
		1 * time.Second,
		2 * time.Second,
		4 * time.Second,
		8 * time.Second,
		10 * time.Second,
		10 * time.Second,
		10 * time.Second,
	}
	for i, expected := range want {
		if got := p.Delay(i + 1); got != expected {
			t.Errorf("Delay(%d) = %v, want %v", i+1, got, expected)
		}
	}
	if got := p.Delay(0); got != time.Second {
		t.Errorf("Delay(0) = %v, want %v", got, time.Second)
	}
	if got := p.Delay(1_000_000); got != 10*time.Second {
		t.Errorf("Delay(1e6) = %v, want saturation at %v", got, 10*time.Second)
	}
}

func TestRetryPolicyJitterBounds(t *testing.T) {
	p := DefaultRetryPolicy()
	p.Rand = fixedRand(0)
	if got := p.Delay(1); got != 800*time.Millisecond {
		t.Errorf("lower jitter bound = %v, want 800ms", got)
	}
	p.Rand = fixedRand(1)
	if got := p.Delay(1); got != 1200*time.Millisecond {
		t.Errorf("upper jitter bound = %v, want 1.2s", got)
	}
}

func TestRetryPolicyClampsToMaxAfterJitter(t *testing.T) {
	p := DefaultRetryPolicy()
	p.Rand = fixedRand(1)
	// 第 5 轮指数项 16s 先截断到 10s，再乘 1.2 后仍需截断到 10s。
	if got := p.Delay(5); got != 10*time.Second {
		t.Fatalf("Delay(5) = %v, want %v", got, 10*time.Second)
	}
}

func TestRetryPolicyZeroValueUsesDefaults(t *testing.T) {
	var p RetryPolicy
	p.Rand = fixedRand(0.5)
	if got, want := p.Delay(1), DefaultRetryBase; got != want {
		t.Fatalf("zero value Delay(1) = %v, want %v", got, want)
	}
}

func TestRetryPolicyCustomShape(t *testing.T) {
	p := RetryPolicy{
		Base:       10 * time.Millisecond,
		Multiplier: 3,
		Max:        100 * time.Millisecond,
		JitterLow:  1,
		JitterHigh: 1,
		Rand:       fixedRand(0),
	}
	want := []time.Duration{10 * time.Millisecond, 30 * time.Millisecond, 90 * time.Millisecond, 100 * time.Millisecond}
	for i, expected := range want {
		if got := p.Delay(i + 1); got != expected {
			t.Errorf("Delay(%d) = %v, want %v", i+1, got, expected)
		}
	}
}
