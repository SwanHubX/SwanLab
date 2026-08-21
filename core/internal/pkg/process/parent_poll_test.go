//go:build !linux && !windows

package process

import (
	"testing"
	"time"
)

func TestWaitForParentExitDetectsInitialMismatch(t *testing.T) {
	getParentPID := func() int { return 2 }
	done := make(chan struct{})
	go func() {
		waitForParentExit(1, make(chan time.Time), getParentPID)
		close(done)
	}()
	select {
	case <-done:
	case <-time.After(time.Second):
		t.Fatal("waitForParentExit did not return on initial mismatch")
	}
}

func TestWaitForParentExitDetectsMismatchOnTick(t *testing.T) {
	const parentPID = 1
	calls := 0
	getParentPID := func() int {
		calls++
		if calls == 1 {
			return parentPID
		}
		return parentPID + 1
	}
	ticks := make(chan time.Time, 1)
	ticks <- time.Time{}

	waitForParentExit(parentPID, ticks, getParentPID)

	if calls != 2 {
		t.Fatalf("getParentPID called %d times, want 2", calls)
	}
}
