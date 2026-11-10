package health

import (
	"context"
	"sync/atomic"
	"testing"
	"time"

	"go.uber.org/goleak"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

// fakeClock 让测试显式推进探测节奏，不用 sleep 制造窗口。
type fakeClock struct {
	created chan chan time.Time
	periods chan time.Duration
}

func newFakeClock() *fakeClock {
	return &fakeClock{
		created: make(chan chan time.Time, 8),
		periods: make(chan time.Duration, 8),
	}
}

func (c *fakeClock) Ticker(d time.Duration) (<-chan time.Time, func()) {
	ch := make(chan time.Time, 8)
	c.created <- ch
	c.periods <- d
	return ch, func() {}
}

// nextPeriod 返回下一个注册的探测周期。
func (c *fakeClock) nextPeriod(t *testing.T) time.Duration {
	t.Helper()
	select {
	case d := <-c.periods:
		return d
	case <-time.After(2 * time.Second):
		t.Fatal("probe ticker was not created")
		return 0
	}
}

// ticks 等待 n 个 prober 注册各自的 ticker。
func (c *fakeClock) ticks(t *testing.T, n int) []chan time.Time {
	t.Helper()
	out := make([]chan time.Time, n)
	for i := range n {
		select {
		case ch := <-c.created:
			out[i] = ch
		case <-time.After(2 * time.Second):
			t.Fatalf("ticker %d/%d not created", i+1, n)
		}
	}
	return out
}

// probe 记录探测调用；seen 无缓冲，prober 必须等测试读取后才继续。
type probe struct {
	ok   atomic.Bool
	seen chan bool
}

func newProbe(ok bool) *probe {
	p := &probe{seen: make(chan bool)}
	p.ok.Store(ok)
	return p
}

func (p *probe) set(ok bool) { p.ok.Store(ok) }

func (p *probe) fn(context.Context) bool {
	ok := p.ok.Load()
	p.seen <- ok
	return ok
}

func (p *probe) wait(t *testing.T) {
	t.Helper()
	select {
	case <-p.seen:
	case <-time.After(2 * time.Second):
		t.Fatal("probe was not called")
	}
}

func tickBoth(t *testing.T, ticks []chan time.Time, control, data *probe) {
	t.Helper()
	for _, ch := range ticks {
		ch <- time.Time{}
	}
	control.wait(t)
	data.wait(t)
}

func expectUpdate(t *testing.T, s *Site) {
	t.Helper()
	select {
	case _, ok := <-s.Updates():
		if !ok {
			t.Fatal("updates channel closed")
		}
	case <-time.After(2 * time.Second):
		t.Fatal("missing state update")
	}
}

func alwaysReachable(context.Context) bool { return true }

func TestProberThresholds(t *testing.T) {
	var reachable atomic.Bool
	var flips []bool
	p := newProber(func(ctx context.Context) bool {
		if _, ok := ctx.Deadline(); !ok {
			t.Error("probe context has no deadline")
		}
		return reachable.Load()
	}, proberOptions{
		timeout:      time.Second,
		degradeAfter: 3,
		recoverAfter: 2,
	}, func(degraded bool) {
		flips = append(flips, degraded)
	})

	for range 2 {
		p.tick(context.Background())
	}
	if p.degraded.Load() {
		t.Fatal("degraded before the failure threshold")
	}
	p.tick(context.Background())
	if !p.degraded.Load() {
		t.Fatal("not degraded at the failure threshold")
	}
	if len(flips) != 1 || !flips[0] {
		t.Fatalf("flips = %v, want [true]", flips)
	}

	reachable.Store(true)
	p.tick(context.Background())
	if !p.degraded.Load() {
		t.Fatal("recovered before the success threshold")
	}
	p.tick(context.Background())
	if p.degraded.Load() {
		t.Fatal("still degraded at the success threshold")
	}
	if len(flips) != 2 || flips[1] {
		t.Fatalf("flips = %v, want [true false]", flips)
	}

	reachable.Store(false)
	for range 3 {
		p.tick(context.Background())
	}
	if !p.degraded.Load() {
		t.Fatal("not degraded after a second outage")
	}
	if len(flips) != 3 || !flips[2] {
		t.Fatalf("flips = %v, want [true false true]", flips)
	}
}

func TestSiteOptionDefaults(t *testing.T) {
	clock := newFakeClock()
	s := newSite(SiteOptions{Key: "defaults", Control: alwaysReachable, Data: alwaysReachable}, clock)
	defer s.stop()

	for range 2 {
		if got := clock.nextPeriod(t); got != DefaultInterval {
			t.Fatalf("default period = %v, want %v", got, DefaultInterval)
		}
	}
	if s.control.timeout != DefaultTimeout {
		t.Fatalf("timeout = %v, want %v", s.control.timeout, DefaultTimeout)
	}
	if s.control.degradeAfter != DefaultDegradeThreshold || s.control.recoverAfter != DefaultRecoverThreshold {
		t.Fatalf("thresholds = %d/%d", s.control.degradeAfter, s.control.recoverAfter)
	}

	s2 := newSite(SiteOptions{
		Key:              "custom",
		Control:          alwaysReachable,
		Data:             alwaysReachable,
		Interval:         time.Minute,
		Timeout:          time.Second,
		DegradeThreshold: 7,
		RecoverThreshold: 4,
	}, clock)
	defer s2.stop()

	for range 2 {
		if got := clock.nextPeriod(t); got != time.Minute {
			t.Fatalf("custom period = %v, want 1m", got)
		}
	}
	if s2.data.timeout != time.Second {
		t.Fatalf("timeout = %v, want 1s", s2.data.timeout)
	}
	if s2.data.degradeAfter != 7 || s2.data.recoverAfter != 4 {
		t.Fatalf("thresholds = %d/%d", s2.data.degradeAfter, s2.data.recoverAfter)
	}
}

func TestSiteAggregatesReasons(t *testing.T) {
	s := newSite(SiteOptions{Key: "site", Control: alwaysReachable, Data: alwaysReachable}, newFakeClock())
	defer s.stop()

	snap := s.Snapshot()
	if snap.State != StateNormal || snap.Reason != "" {
		t.Fatalf("initial snapshot = %+v", snap)
	}

	s.setReason(ReasonProbeData, true)
	expectUpdate(t, s)
	if snap := s.Snapshot(); snap.State != StateDegraded || snap.Reason != ReasonProbeData {
		t.Fatalf("data degraded snapshot = %+v", snap)
	}

	s.setReason(ReasonProbeControl, true)
	expectUpdate(t, s)
	if snap := s.Snapshot(); snap.State != StateDegraded || snap.Reason != ReasonProbeControl {
		t.Fatalf("control should take precedence, snapshot = %+v", snap)
	}

	s.setReason(ReasonProbeControl, false)
	expectUpdate(t, s)
	if snap := s.Snapshot(); snap.Reason != ReasonProbeData {
		t.Fatalf("snapshot = %+v, want data reason", snap)
	}

	s.setReason(ReasonProbeData, false)
	expectUpdate(t, s)
	if snap := s.Snapshot(); snap.State != StateNormal || snap.Reason != "" {
		t.Fatalf("recovered snapshot = %+v", snap)
	}

	s.Report(ReasonWrite, true)
	expectUpdate(t, s)
	if snap := s.Snapshot(); snap.State != StateDegraded || snap.Reason != ReasonWrite {
		t.Fatalf("write degraded snapshot = %+v", snap)
	}

	s.Report(ReasonWrite, false)
	expectUpdate(t, s)
	if snap := s.Snapshot(); snap.State != StateNormal {
		t.Fatalf("write recovered snapshot = %+v", snap)
	}

	// 外部只能注入写入侧原因，否则会与 prober 的状态脱节。
	func() {
		defer func() {
			if recover() == nil {
				t.Fatal("Report(ReasonProbeData) should panic")
			}
		}()
		s.Report(ReasonProbeData, true)
	}()
}

func TestRegistrySharesSites(t *testing.T) {
	clock := newFakeClock()
	reg := NewRegistry(clock)
	control, data := newProbe(true), newProbe(true)

	a := reg.Acquire(SiteOptions{Key: "site", Control: control.fn, Data: data.fn})
	ticks := clock.ticks(t, 2)

	b := reg.Acquire(SiteOptions{Key: "site", Control: alwaysReachable, Data: alwaysReachable})
	if a != b {
		t.Fatal("same key should share one site")
	}
	if got := reg.Len(); got != 1 {
		t.Fatalf("sites = %d, want 1", got)
	}

	reg.Release(a)
	if got := reg.Len(); got != 1 {
		t.Fatalf("sites = %d after one release, want 1", got)
	}
	// 仍有引用的站点必须继续探测。
	tickBoth(t, ticks, control, data)

	reg.Release(b)
	if got := reg.Len(); got != 0 {
		t.Fatalf("sites = %d after all releases, want 0", got)
	}
	if _, ok := <-a.Updates(); ok {
		t.Fatal("updates channel should be closed after release")
	}

	// 释放在前一次的 ID 上重复调用不得影响新实例。
	reg.Release(a)
	c := reg.Acquire(SiteOptions{Key: "site", Control: alwaysReachable, Data: alwaysReachable})
	if c == a {
		t.Fatal("released site must not be reused")
	}
	reg.Release(c)
	reg.Release(nil)
}

func TestRegistryMergesLaterInitialDegraded(t *testing.T) {
	clock := newFakeClock()
	reg := NewRegistry(clock)
	control, data := newProbe(true), newProbe(true)

	a := reg.Acquire(SiteOptions{Key: "site", Control: control.fn, Data: data.fn})
	clock.ticks(t, 2)
	b := reg.Acquire(SiteOptions{
		Key:                 "site",
		Control:             alwaysReachable,
		Data:                alwaysReachable,
		InitialDataDegraded: true,
	})
	if a != b {
		t.Fatal("same key should share one site")
	}
	expectUpdate(t, a)
	snap := a.Snapshot()
	if snap.State != StateDegraded || snap.Reason != ReasonProbeData {
		t.Fatalf("snapshot = %+v, want degraded by later data precheck", snap)
	}
	reg.Release(a)
	reg.Release(b)
}

func TestRegistryProbeDrivesState(t *testing.T) {
	clock := newFakeClock()
	control, data := newProbe(true), newProbe(false)

	reg := NewRegistry(clock)
	site := reg.Acquire(SiteOptions{
		Key:              "site",
		Control:          control.fn,
		Data:             data.fn,
		DegradeThreshold: 3,
		RecoverThreshold: 2,
	})
	defer reg.Release(site)
	ticks := clock.ticks(t, 2)

	for range 3 {
		tickBoth(t, ticks, control, data)
	}
	expectUpdate(t, site)
	snap := site.Snapshot()
	if snap.State != StateDegraded || snap.Reason != ReasonProbeData {
		t.Fatalf("snapshot = %+v, want degraded by data", snap)
	}

	data.set(true)
	for range 2 {
		tickBoth(t, ticks, control, data)
	}
	expectUpdate(t, site)
	snap = site.Snapshot()
	if snap.State != StateNormal || snap.Reason != "" {
		t.Fatalf("snapshot = %+v, want normal", snap)
	}
}

func TestAcquireRejectsInvalidOptions(t *testing.T) {
	reg := NewRegistry(newFakeClock())
	cases := map[string]SiteOptions{
		"missing key":   {Control: alwaysReachable, Data: alwaysReachable},
		"missing probe": {Key: "site", Data: alwaysReachable},
	}
	for name, opts := range cases {
		t.Run(name, func(t *testing.T) {
			defer func() {
				if recover() == nil {
					t.Fatal("expected panic")
				}
			}()
			reg.Acquire(opts)
		})
	}
	if got := reg.Len(); got != 0 {
		t.Fatalf("sites = %d, want 0", got)
	}
}

func TestSiteReportAfterStop(t *testing.T) {
	reg := NewRegistry(newFakeClock())
	site := reg.Acquire(SiteOptions{Key: "site", Control: alwaysReachable, Data: alwaysReachable})
	reg.Release(site)

	site.Report(ReasonWrite, true)
	site.Report(ReasonWrite, false)
	if snap := site.Snapshot(); snap.State != StateNormal {
		t.Fatalf("snapshot after stop = %+v", snap)
	}
	reg.Release(site)
}

func TestSiteInitialDegradedRecovers(t *testing.T) {
	clock := newFakeClock()
	control, data := newProbe(true), newProbe(true)
	s := newSite(SiteOptions{
		Key:                 "seed",
		Control:             control.fn,
		Data:                data.fn,
		InitialDataDegraded: true,
	}, clock)
	defer s.stop()
	ticks := clock.ticks(t, 2)

	snap := s.Snapshot()
	if snap.State != StateDegraded || snap.Reason != ReasonProbeData {
		t.Fatalf("initial snapshot = %+v, want degraded by data", snap)
	}
	// 消费构造期种下的通知，后面的断言只等待恢复事件。
	expectUpdate(t, s)

	for range 2 {
		tickBoth(t, ticks, control, data)
	}
	expectUpdate(t, s)
	if snap := s.Snapshot(); snap.State != StateNormal || snap.Reason != "" {
		t.Fatalf("snapshot after recovery = %+v, want normal", snap)
	}
}
