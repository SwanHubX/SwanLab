package heartbeat

import (
	"context"
	"errors"
	"testing"
	"time"

	"go.uber.org/goleak"
)

const heartbeatWait = time.Second

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

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
	hb := New("exp-1", func(_ context.Context, id string) error {
		calls <- id
		return nil
	})
	hb.Start(context.Background())
	t.Cleanup(func() { hb.Stop(context.Background()) })

	if got := waitBeat(t, calls); got != "exp-1" {
		t.Fatalf("heartbeat experiment = %q, want exp-1", got)
	}
}

func TestHeartbeatRepeats(t *testing.T) {
	calls := make(chan string, 4)
	hb := New("exp-1", func(_ context.Context, id string) error {
		calls <- id
		return nil
	})
	hb.interval = 5 * time.Millisecond
	hb.Start(context.Background())
	t.Cleanup(func() { hb.Stop(context.Background()) })

	waitBeat(t, calls)
	waitBeat(t, calls)
}

func TestHeartbeatStopIsIdempotent(_ *testing.T) {
	hb := New("exp-1", func(context.Context, string) error { return nil })
	hb.Start(context.Background())
	hb.Stop(context.Background())
	hb.Stop(context.Background())
}

func TestHeartbeatStopBeforeStart(t *testing.T) {
	hb := New("exp-1", func(context.Context, string) error {
		t.Error("heartbeat must not send after stop before start")
		return nil
	})
	hb.Stop(context.Background())
	hb.Start(context.Background())
	hb.Stop(context.Background())
}

func TestHeartbeatCanceledContextDoesNotSend(t *testing.T) {
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	hb := New("exp-1", func(context.Context, string) error {
		t.Error("heartbeat must not send on a canceled context")
		return nil
	})
	hb.Start(ctx)
	hb.Stop(context.Background())
}

func TestHeartbeatSendErrorIsBestEffort(t *testing.T) {
	called := make(chan struct{}, 1)
	hb := New("exp-1", func(context.Context, string) error {
		called <- struct{}{}
		return errors.New("cloud unavailable")
	})
	hb.Start(context.Background())
	t.Cleanup(func() { hb.Stop(context.Background()) })

	select {
	case <-called:
	case <-time.After(heartbeatWait):
		t.Fatal("timed out waiting for heartbeat")
	}
}

// Stop 必须先取消在途请求，否则慢心跳会把调用方的等待拖到 HTTP 超时。
func TestHeartbeatStopCancelsInFlightSend(t *testing.T) {
	entered := make(chan struct{})
	hb := New("exp-1", func(ctx context.Context, _ string) error {
		close(entered)
		<-ctx.Done()
		return ctx.Err()
	})
	hb.Start(context.Background())

	select {
	case <-entered:
	case <-time.After(heartbeatWait):
		t.Fatal("heartbeat did not start")
	}

	stopped := make(chan struct{})
	go func() {
		hb.Stop(context.Background())
		close(stopped)
	}()
	select {
	case <-stopped:
	case <-time.After(heartbeatWait):
		t.Fatal("stop did not cancel the in-flight heartbeat")
	}
}
