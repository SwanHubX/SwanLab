package manager

import (
	"context"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/cron/health"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

func okProbe(context.Context) bool { return true }

func TestWriteFailureDegradesSite(t *testing.T) {
	reg := health.NewRegistry(nil)
	s := newSession(sessionOptions{})
	s.healthReg = reg
	s.site = reg.Acquire(health.SiteOptions{Key: "unit", Server: okProbe, House: okProbe})
	t.Cleanup(s.closeHealth)

	s.onWriteFailure(2)
	if snap := s.site.Snapshot(); snap.State != health.StateNormal {
		t.Fatalf("snapshot = %+v, want normal before threshold", snap)
	}

	s.onWriteFailure(3)
	snap := s.site.Snapshot()
	if snap.State != health.StateDegraded || snap.Reason != health.ReasonWrite {
		t.Fatalf("snapshot = %+v, want degraded by write", snap)
	}

	s.onWriteSuccess()
	if snap := s.site.Snapshot(); snap.State != health.StateNormal {
		t.Fatalf("snapshot = %+v, want normal after success", snap)
	}
}

// confirmAll 立即确认整批记录，用于构造空闲 uploader。
type confirmAll struct{}

func (confirmAll) Upload(_ context.Context, records []*recordv1.Record) uploader.Result {
	return uploader.Result{Confirmed: records}
}

func newIdleUploader(t *testing.T) *uploader.Uploader {
	t.Helper()
	up, err := uploader.New(uploader.Options{
		In:       make(chan *recordv1.Record),
		Budget:   budget.New(budget.Limits{}),
		Senders:  map[record.Type]uploader.Sender{record.TypeScalar: confirmAll{}},
		Interval: time.Hour,
	})
	if err != nil {
		t.Fatalf("uploader.New: %v", err)
	}
	t.Cleanup(func() { _ = up.Close() })
	return up
}

func TestReconcileHealthSuspendsOnlyWithLocalCopy(t *testing.T) {
	reg := health.NewRegistry(nil)
	site := reg.Acquire(health.SiteOptions{
		Key:                  "unit",
		Server:               okProbe,
		House:                okProbe,
		InitialHouseDegraded: true,
	})
	t.Cleanup(func() { reg.Release(site) })

	up := newIdleUploader(t)
	s := &Session{cfg: config.Defaults(), uploader: up, site: site}

	s.reconcileHealth()
	if !up.Suspended() {
		t.Fatal("degraded site with a local copy should suspend uploads")
	}

	s.storeDegraded.Store(true)
	s.reconcileHealth()
	if up.Suspended() {
		t.Fatal("store degraded must not keep uploads suspended")
	}

	s.storeDegraded.Store(false)
	s.reconcileHealth()
	if !up.Suspended() {
		t.Fatal("healthy store should suspend again while the site is degraded")
	}

	s.cfg.SkipStore = true
	s.reconcileHealth()
	if up.Suspended() {
		t.Fatal("skip_store must not suspend uploads")
	}
}

func TestWriteDegradeDoesNotSuspend(t *testing.T) {
	reg := health.NewRegistry(nil)
	site := reg.Acquire(health.SiteOptions{Key: "unit", Server: okProbe, House: okProbe})
	t.Cleanup(func() { reg.Release(site) })

	up := newIdleUploader(t)
	s := &Session{cfg: config.Defaults(), uploader: up, site: site}

	site.Report(health.ReasonWrite, true)
	s.reconcileHealth()
	if up.Suspended() {
		t.Fatal("write degradation must keep uploads running so recovery can happen")
	}
}

func TestSharedSiteHealthReconcilesAllSessions(t *testing.T) {
	reg := health.NewRegistry(nil)
	houseRecovery := make(chan struct{})
	houseProbe := func(context.Context) bool {
		select {
		case <-houseRecovery:
			return true
		default:
			return false
		}
	}
	opts := health.SiteOptions{
		Key:              "shared",
		Interval:         2 * time.Millisecond,
		DegradeThreshold: 1,
		RecoverThreshold: 1,
		Server:           okProbe,
		House:            houseProbe,
	}
	a := newSession(sessionOptions{})
	a.cfg = config.Defaults()
	a.uploader = newIdleUploader(t)
	b := newSession(sessionOptions{})
	b.cfg = config.Defaults()
	b.uploader = newIdleUploader(t)
	t.Cleanup(func() {
		_ = a.uploader.Close()
		_ = b.uploader.Close()
		a.closeHealth()
		b.closeHealth()
	})
	a.attachHealth(reg, opts)
	b.attachHealth(reg, opts)

	wait := func(want bool) {
		t.Helper()
		ticker := time.NewTicker(time.Millisecond)
		defer ticker.Stop()
		timeout := time.NewTimer(2 * time.Second)
		defer timeout.Stop()
		for {
			if a.uploader.Suspended() == want && b.uploader.Suspended() == want {
				return
			}
			select {
			case <-ticker.C:
			case <-timeout.C:
				t.Fatalf("upload suspended states = %t/%t, want %t",
					a.uploader.Suspended(), b.uploader.Suspended(), want)
			}
		}
	}
	wait(true)
	close(houseRecovery)
	wait(false)
}
