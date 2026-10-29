package uploader

import (
	"sync"
	"testing"
	"time"
)

// throttleRecorder 收集节流器输出，替代真实终端写入。
type throttleRecorder struct {
	mu    sync.Mutex
	warns []string
	infos []string
}

func (r *throttleRecorder) warn(msg string) {
	r.mu.Lock()
	defer r.mu.Unlock()
	r.warns = append(r.warns, msg)
}

func (r *throttleRecorder) info(msg string) {
	r.mu.Lock()
	defer r.mu.Unlock()
	r.infos = append(r.infos, msg)
}

func (r *throttleRecorder) snapshot() (warns, infos []string) {
	r.mu.Lock()
	defer r.mu.Unlock()
	return append([]string(nil), r.warns...), append([]string(nil), r.infos...)
}

func newThrottleFixture() (*WarningThrottle, *throttleRecorder, *fakeClock) {
	clock := newFakeClock()
	recorder := &throttleRecorder{}
	throttle := newWarningThrottle(clock.now, DefaultWarningInterval)
	throttle.warnf = recorder.warn
	throttle.infof = recorder.info
	return throttle, recorder, clock
}

func TestWarningThrottleFirstWarnIsImmediate(t *testing.T) {
	throttle, recorder, clock := newThrottleFixture()

	throttle.Warn()
	warns, _ := recorder.snapshot()
	if len(warns) != 1 || warns[0] != uploadFailureWarning {
		t.Fatalf("warns = %v, want one failure warning", warns)
	}

	clock.advance(DefaultWarningInterval - time.Second)
	throttle.Warn()
	if warns, _ = recorder.snapshot(); len(warns) != 1 {
		t.Fatalf("warns within interval = %d, want 1", len(warns))
	}

	clock.advance(time.Second)
	throttle.Warn()
	if warns, _ = recorder.snapshot(); len(warns) != 2 {
		t.Fatalf("warns at interval boundary = %d, want 2", len(warns))
	}
}

func TestWarningThrottleReset(t *testing.T) {
	throttle, recorder, clock := newThrottleFixture()

	// 未发生故障时 Reset 不提示。
	throttle.Reset()
	if _, infos := recorder.snapshot(); len(infos) != 0 {
		t.Fatalf("reset without failure should stay silent: %v", infos)
	}

	throttle.Warn()
	clock.advance(time.Second) // 处于故障期，但未到下次告警

	throttle.Reset()
	_, infos := recorder.snapshot()
	if len(infos) != 1 || infos[0] != uploadRecoveredNotice {
		t.Fatalf("infos = %v, want one recovery notice", infos)
	}

	throttle.Reset()
	if _, infos = recorder.snapshot(); len(infos) != 1 {
		t.Fatalf("second reset should stay silent: %v", infos)
	}

	// Reset 复位告警窗口，下一次失败立即告警。
	throttle.Warn()
	if warns, _ := recorder.snapshot(); len(warns) != 2 {
		t.Fatalf("warn after reset = %d, want 2", len(warns))
	}
}

func TestWarningThrottleConcurrent(t *testing.T) {
	throttle, recorder, _ := newThrottleFixture()

	var wg sync.WaitGroup
	for i := 0; i < 50; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			throttle.Warn()
		}()
	}
	wg.Wait()

	if warns, _ := recorder.snapshot(); len(warns) != 1 {
		t.Fatalf("concurrent warns = %d, want 1", len(warns))
	}
}
