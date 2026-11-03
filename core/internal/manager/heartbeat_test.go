package manager

import (
	"context"
	"errors"
	"testing"
	"time"
)

const heartbeatWait = time.Second

func waitBeat(t *testing.T, calls <-chan string) string {
	t.Helper()
	select {
	case id := <-calls:
		return id
	case <-time.After(heartbeatWait):
		t.Fatal("timed out waiting for heartbeat")
		return ""
	}
}

func TestHeartbeatSendsImmediately(t *testing.T) {
	calls := make(chan string, 1)
	hb := newHeartbeat("exp-1", func(_ context.Context, id string) error {
		calls <- id
		return nil
	})
	hb.start(context.Background())
	t.Cleanup(hb.stop)

	if got := waitBeat(t, calls); got != "exp-1" {
		t.Fatalf("heartbeat experiment = %q, want exp-1", got)
	}
}

func TestHeartbeatRepeats(t *testing.T) {
	calls := make(chan string, 4)
	hb := newHeartbeat("exp-1", func(_ context.Context, id string) error {
		calls <- id
		return nil
	})
	hb.interval = 5 * time.Millisecond
	hb.start(context.Background())
	t.Cleanup(hb.stop)

	waitBeat(t, calls)
	waitBeat(t, calls)
}

func TestHeartbeatStopIsIdempotent(_ *testing.T) {
	hb := newHeartbeat("exp-1", func(context.Context, string) error { return nil })
	hb.start(context.Background())
	hb.stop()
	hb.stop()
}

func TestHeartbeatStopBeforeStart(t *testing.T) {
	hb := newHeartbeat("exp-1", func(context.Context, string) error {
		t.Error("heartbeat must not send after stop before start")
		return nil
	})
	hb.stop()
	hb.start(context.Background())
	hb.stop()
}

func TestHeartbeatCanceledContextDoesNotSend(t *testing.T) {
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	hb := newHeartbeat("exp-1", func(context.Context, string) error {
		t.Error("heartbeat must not send on a canceled context")
		return nil
	})
	hb.start(ctx)
	hb.stop()
}

func TestHeartbeatSendErrorIsBestEffort(t *testing.T) {
	called := make(chan struct{}, 1)
	hb := newHeartbeat("exp-1", func(context.Context, string) error {
		called <- struct{}{}
		return errors.New("cloud unavailable")
	})
	hb.start(context.Background())
	t.Cleanup(hb.stop)

	select {
	case <-called:
	case <-time.After(heartbeatWait):
		t.Fatal("timed out waiting for heartbeat")
	}
}
