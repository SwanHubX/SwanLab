package manager

import (
	"context"
	"encoding/json"
	"errors"
	"net/http"
	"path/filepath"
	"strings"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/manager/probe"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func fastAttachTiming() attachTiming {
	return attachTiming{
		initialBackoff: 5 * time.Millisecond,
		maxBackoff:     20 * time.Millisecond,
		jitterMax:      time.Millisecond,
		attemptTimeout: 2 * time.Second,
	}
}

func newInstantTestManager() *Manager {
	mgr := New(Options{
		HTTP: client.Options{
			Timeout: time.Second,
			Sleep:   func(context.Context, time.Duration) error { return nil },
		},
	})
	mgr.attachTiming = fastAttachTiming()
	return mgr
}

func waitAttachDone(t *testing.T, s *Session) {
	t.Helper()
	timer := time.NewTimer(5 * time.Second)
	defer timer.Stop()
	select {
	case <-s.attachDone:
	case <-timer.C:
		state, path, err := s.cloudStatus()
		t.Fatalf("timed out waiting for attachDone (state=%v path=%q err=%v)", state, path, err)
	}
}

func waitManagerCond(t *testing.T, desc string, cond func() bool) {
	t.Helper()
	if cond() {
		return
	}
	ticker := time.NewTicker(2 * time.Millisecond)
	defer ticker.Stop()
	timeout := time.NewTimer(5 * time.Second)
	defer timeout.Stop()
	for {
		select {
		case <-ticker.C:
			if cond() {
				return
			}
		case <-timeout.C:
			t.Fatalf("timed out waiting for %s", desc)
		}
	}
}

func serveMetricsUploadOK(cloud *clienttest.Cloud) <-chan []byte {
	cloud.ServeJSON(http.MethodPost, "/api/projects/alice/demo/series", http.StatusOK, map[string]any{"ok": true})
	ch := make(chan []byte, 32)
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, r *http.Request) {
		var payload map[string]any
		_ = json.NewDecoder(r.Body).Decode(&payload)
		raw, _ := json.Marshal(payload)
		select {
		case ch <- raw:
		default:
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	return ch
}

func newDirectAttachSession(
	t *testing.T,
	mgr *Manager,
	cloud *clienttest.Cloud,
	serverDown bool,
	houseDown bool,
) (*Session, *runv1.StartRecord, func()) {
	t.Helper()
	cfg := startSettings(t, cloud)
	cfg.RecoverOnline = true
	cfg.HealthProbeInterval = 0.005
	cfg.HealthProbeTimeout = 0.2
	cfg.HealthDegradeThreshold = 1
	cfg.HealthRecoverThreshold = 1

	run, err := runctx.New(cfg)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	serverProbe := probe.Site(cfg.APIHost, probe.ServerPath)
	houseProbe := probe.Site(cfg.APIHost, probe.HousePath)
	healthOpts := buildHealthOpts(cfg, serverProbe, houseProbe, serverDown, houseDown)

	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.run = run
	s.metrics = metrics.New()
	s.counter = counter.New(0)
	s.epoch = counter.New(0)
	s.initDeferredAttach()

	if err := s.open(sessionInit{
		cfg:        cfg,
		run:        run,
		healthReg:  mgr.healthReg,
		healthOpts: healthOpts,
	}); err != nil {
		s.rollback()
		t.Fatalf("s.open: %v", err)
	}
	start := &runv1.StartRecord{
		Id:        cfg.RunID,
		Workspace: "alice",
		Project:   "demo",
		Name:      "exp-recover",
		Color:     "#abcdef",
	}
	if err := s.commitStart(context.Background(), start); err != nil {
		s.rollback()
		t.Fatalf("commitStart: %v", err)
	}
	cleanup := func() {
		ctx, cancel := context.WithTimeout(context.Background(), 2*time.Second)
		defer cancel()
		_ = s.shutdown(ctx)
		if up := s.uploader(); up != nil {
			_ = up.Close()
		}
	}
	return s, start, cleanup
}

func TestAttachRetriesOn5xxAndAttachesOnce(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, probe.ServerPath, http.StatusOK, map[string]any{"ok": true})
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
	serveAuth(cloud)
	serveProject(cloud, "alice", "demo", 1)
	serveHeartbeat(cloud, "exp-1")
	uploaded := serveMetricsUploadOK(cloud)

	var expCalls atomic.Int32
	cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		n := expCalls.Add(1)
		if n <= 2 {
			clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]any{"message": "busy"})
			return
		}
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		clienttest.WriteJSON(w, http.StatusCreated, map[string]any{
			"cuid": "exp-1",
			"slug": "slug-1",
			"name": body["name"],
		})
	})

	mgr := newInstantTestManager()
	s, start, cleanup := newDirectAttachSession(t, mgr, cloud, false, false)
	defer cleanup()

	// 接入前写入一条标量，应当仅落盘而不上传。
	id, err := mgr.install(s)
	if err != nil {
		t.Fatalf("install: %v", err)
	}
	if err := mgr.UpsertScalars(context.Background(), id, []*datav1.ScalarRecord{
		{
			Key:       "loss",
			Step:      1,
			Type:      columnv1.ColumnType_COLUMN_TYPE_SCALAR,
			Value:     &datav1.ScalarValue{Number: 0.25},
			Timestamp: timestamppb.New(time.Unix(1_700_000_000, 0)),
		},
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}

	s.startAttach(attachParams{
		factory:    mgr.factory,
		start:      start,
		resumeMode: api.ResumeNever,
		timing:     fastAttachTiming(),
	})

	waitAttachDone(t, s)
	state, path, cloudErr := s.cloudStatus()
	if state != operationv1.CloudState_CLOUD_STATE_ATTACHED || path != "/alice/demo/slug-1" || cloudErr != nil {
		t.Fatalf("cloudStatus = (%v, %q, %v), want ATTACHED /alice/demo/slug-1", state, path, cloudErr)
	}
	if got := expCalls.Load(); got != 3 {
		t.Fatalf("CreateOrResumeExperiment calls = %d, want 3 (2 failures + 1 created)", got)
	}
	params, ok := s.run.Online()
	if !ok || params.ExperimentID != "exp-1" || params.ProjectID != "proj-1" || params.Username != "alice" {
		t.Fatalf("run.Online() = (%+v, %v), want attached identity", params, ok)
	}

	// 验证历史记录通过 backfill 补传到云端。
	timer := time.NewTimer(5 * time.Second)
	defer timer.Stop()
	select {
	case <-uploaded:
	case <-timer.C:
		t.Fatal("timed out waiting for buffered scalar upload after attach")
	}

	snap := s.snapshot()
	if snap.GetCloud().GetState() != operationv1.CloudState_CLOUD_STATE_ATTACHED || snap.GetCloud().GetPath() != "/alice/demo/slug-1" {
		t.Fatalf("snapshot cloud = %+v", snap.GetCloud())
	}
}

func TestAttachWakesOnServerRecovery(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var serverUp atomic.Bool
	cloud.Handle(http.MethodGet, probe.ServerPath, func(w http.ResponseWriter, _ *http.Request) {
		if !serverUp.Load() {
			clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]any{"ok": false})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
	serveAuth(cloud)
	serveProject(cloud, "alice", "demo", 1)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "slug-wake")
	serveHeartbeat(cloud, "exp-1")
	uploaded := serveMetricsUploadOK(cloud)

	mgr := newInstantTestManager()
	s, start, cleanup := newDirectAttachSession(t, mgr, cloud, true, false)
	defer cleanup()

	id, err := mgr.install(s)
	if err != nil {
		t.Fatalf("install: %v", err)
	}
	if err := mgr.UpsertScalars(context.Background(), id, []*datav1.ScalarRecord{
		{
			Key:       "acc",
			Step:      1,
			Type:      columnv1.ColumnType_COLUMN_TYPE_SCALAR,
			Value:     &datav1.ScalarValue{Number: 0.9},
			Timestamp: timestamppb.New(time.Unix(1_700_000_000, 0)),
		},
	}); err != nil {
		t.Fatalf("UpsertScalars: %v", err)
	}

	// 长退避期间由站点恢复通知立即唤醒接入，而不必等长退避到期。
	s.startAttach(attachParams{
		factory:    mgr.factory,
		start:      start,
		resumeMode: api.ResumeNever,
		timing: attachTiming{
			initialBackoff: 10 * time.Second,
			maxBackoff:     30 * time.Second,
			jitterMax:      time.Millisecond,
			attemptTimeout: 2 * time.Second,
		},
	})

	if snap := s.snapshot(); snap.GetCloud().GetState() != operationv1.CloudState_CLOUD_STATE_PENDING {
		t.Fatalf("initial cloud state = %v, want PENDING", snap.GetCloud().GetState())
	}
	if cloud.CountFor(http.MethodPost, "/api/project/alice/demo/experiment") != 0 {
		t.Fatal("must not call experiment API while server probe is degraded")
	}

	serverUp.Store(true)
	waitAttachDone(t, s)

	state, path, cloudErr := s.cloudStatus()
	if state != operationv1.CloudState_CLOUD_STATE_ATTACHED || path != "/alice/demo/slug-wake" || cloudErr != nil {
		t.Fatalf("cloudStatus = (%v, %q, %v)", state, path, cloudErr)
	}
	timer := time.NewTimer(5 * time.Second)
	defer timer.Stop()
	select {
	case <-uploaded:
	case <-timer.C:
		t.Fatal("timed out waiting for upload after server recovery")
	}
}

func TestAttachWhileHouseStillDegradedKeepsUploaderSuspended(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, probe.ServerPath, http.StatusOK, map[string]any{"ok": true})
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusServiceUnavailable, map[string]any{"ok": false})
	serveAuth(cloud)
	serveProject(cloud, "alice", "demo", 1)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "slug-house-down")
	serveHeartbeat(cloud, "exp-1")

	mgr := newInstantTestManager()
	s, start, cleanup := newDirectAttachSession(t, mgr, cloud, false, true)
	defer cleanup()

	s.startAttach(attachParams{
		factory:    mgr.factory,
		start:      start,
		resumeMode: api.ResumeNever,
		timing:     fastAttachTiming(),
	})

	waitAttachDone(t, s)
	state, path, cloudErr := s.cloudStatus()
	if state != operationv1.CloudState_CLOUD_STATE_ATTACHED || path != "/alice/demo/slug-house-down" || cloudErr != nil {
		t.Fatalf("cloudStatus = (%v, %q, %v)", state, path, cloudErr)
	}
	if up := s.uploader(); up == nil || !up.Suspended() {
		t.Fatalf("uploader suspended = %v, want true while house probe is degraded", up != nil && up.Suspended())
	}
}

func TestAttachAuthFailurePermanentAfterThreeConsecutive(t *testing.T) {
	t.Run("non-consecutive 401 resets counter and attaches", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		cloud.ServeJSON(http.MethodGet, probe.ServerPath, http.StatusOK, map[string]any{"ok": true})
		cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
		serveHeartbeat(cloud, "exp-1")

		var verifyCalls atomic.Int32
		cloud.Handle(http.MethodGet, "/api/auth/verify", func(w http.ResponseWriter, _ *http.Request) {
			switch verifyCalls.Add(1) {
			case 1, 2:
				clienttest.WriteJSON(w, http.StatusUnauthorized, map[string]any{"code": "unauth", "message": "temporary 401"})
			case 3:
				clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]any{"code": "busy", "message": "server busy"})
			case 4, 5:
				clienttest.WriteJSON(w, http.StatusUnauthorized, map[string]any{"code": "unauth", "message": "temporary 401"})
			default:
				clienttest.WriteJSON(w, http.StatusOK, map[string]any{"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"})
			}
		})
		serveProject(cloud, "alice", "demo", 1)
		serveNewExperiment(cloud, "alice", "demo", "exp-1", "slug-reset")

		mgr := newInstantTestManager()
		s, start, cleanup := newDirectAttachSession(t, mgr, cloud, false, false)
		defer cleanup()

		s.startAttach(attachParams{
			factory:    mgr.factory,
			start:      start,
			resumeMode: api.ResumeNever,
			timing:     fastAttachTiming(),
		})

		waitAttachDone(t, s)
		state, _, cloudErr := s.cloudStatus()
		if state != operationv1.CloudState_CLOUD_STATE_ATTACHED || cloudErr != nil {
			t.Fatalf("cloudStatus = (%v, %v), want ATTACHED", state, cloudErr)
		}
	})

	for _, statusCode := range []int{http.StatusUnauthorized, http.StatusForbidden} {
		t.Run(http.StatusText(statusCode), func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			cloud.ServeJSON(http.MethodGet, probe.ServerPath, http.StatusOK, map[string]any{"ok": true})
			cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
			cloud.ServeJSON(http.MethodGet, "/api/auth/verify", statusCode, map[string]any{
				"code":    "auth_error",
				"message": "invalid api key detail",
			})

			mgr := newInstantTestManager()
			s, start, cleanup := newDirectAttachSession(t, mgr, cloud, false, false)
			defer cleanup()

			s.startAttach(attachParams{
				factory:    mgr.factory,
				start:      start,
				resumeMode: api.ResumeNever,
				timing:     fastAttachTiming(),
			})

			waitAttachDone(t, s)
			state, path, cloudErr := s.cloudStatus()
			if state != operationv1.CloudState_CLOUD_STATE_FAILED || path != "" || cloudErr == nil {
				t.Fatalf("cloudStatus = (%v, %q, %v), want FAILED", state, path, cloudErr)
			}
			if !strings.Contains(cloudErr.Error(), "invalid api key detail") {
				t.Fatalf("cloudErr = %v, want server reason preserved", cloudErr)
			}
			if got := cloud.CountFor(http.MethodGet, "/api/auth/verify"); got != attachAuthFailLimit {
				t.Fatalf("auth verify calls = %d, want %d", got, attachAuthFailLimit)
			}
			if snap := s.snapshot(); !strings.Contains(snap.GetCloud().GetError(), "invalid api key detail") {
				t.Fatalf("snapshot cloud = %+v, want error detail", snap.GetCloud())
			}
		})
	}
}

func TestAttachImmediatePermanentFailure(t *testing.T) {
	cases := []struct {
		name     string
		resume   api.ResumeMode
		setupAPI func(*clienttest.Cloud)
		wantErr  string
		wantIs   error
		check    func(*testing.T, *clienttest.Cloud)
	}{
		{
			name:   "bad request on create experiment",
			resume: api.ResumeNever,
			setupAPI: func(c *clienttest.Cloud) {
				serveAuth(c)
				serveProject(c, "alice", "demo", 1)
				c.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusBadRequest, map[string]any{
					"code":    "bad_request",
					"message": "experiment name invalid",
				})
			},
			wantErr: "experiment name invalid",
		},
		{
			name:   "resume must missing cloud experiment",
			resume: api.ResumeMust,
			setupAPI: func(c *clienttest.Cloud) {
				serveAuth(c)
				serveProject(c, "alice", "demo", 1)
				c.ServeJSON(http.MethodGet, "/api/project/alice/demo/runs/run-1", http.StatusNotFound, map[string]any{
					"code":    "not_found",
					"message": "experiment does not exist",
				})
			},
			wantErr: "does not exist",
			check: func(t *testing.T, c *clienttest.Cloud) {
				t.Helper()
				if got := c.CountFor(http.MethodGet, "/api/project/alice/demo/runs/run-1"); got != 1 {
					t.Fatalf("resume must lookup count = %d, want 1", got)
				}
			},
		},
		{
			name:   "invalid response incomplete project metadata",
			resume: api.ResumeNever,
			setupAPI: func(c *clienttest.Cloud) {
				serveAuth(c)
				c.ServeJSON(http.MethodGet, "/api/project/alice/demo", http.StatusOK, map[string]any{
					"name": "demo",
				})
			},
			wantErr: "project response missing cuid",
			wantIs:  client.ErrInvalidResponse,
			check: func(t *testing.T, c *clienttest.Cloud) {
				t.Helper()
				if got := c.CountFor(http.MethodGet, "/api/project/alice/demo"); got != 1 {
					t.Fatalf("project lookup count = %d, want 1", got)
				}
			},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			cloud.ServeJSON(http.MethodGet, probe.ServerPath, http.StatusOK, map[string]any{"ok": true})
			cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
			tc.setupAPI(cloud)

			mgr := newInstantTestManager()
			s, start, cleanup := newDirectAttachSession(t, mgr, cloud, false, false)
			defer cleanup()

			s.startAttach(attachParams{
				factory:    mgr.factory,
				start:      start,
				resumeMode: tc.resume,
				timing:     fastAttachTiming(),
			})

			waitAttachDone(t, s)
			state, _, cloudErr := s.cloudStatus()
			if state != operationv1.CloudState_CLOUD_STATE_FAILED || cloudErr == nil {
				t.Fatalf("cloudStatus = (%v, %v), want FAILED", state, cloudErr)
			}
			if tc.wantIs != nil && !errors.Is(cloudErr, tc.wantIs) {
				t.Fatalf("cloudErr = %v, want errors.Is %v", cloudErr, tc.wantIs)
			}
			if !strings.Contains(cloudErr.Error(), tc.wantErr) {
				t.Fatalf("cloudErr = %v, want containing %q", cloudErr, tc.wantErr)
			}
			if tc.check != nil {
				tc.check(t, cloud)
			}
		})
	}
}

func TestAttachDoesNotStartBackfillAfterDraining(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, probe.ServerPath, http.StatusOK, map[string]any{"ok": true})
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
	serveAuth(cloud)
	serveProject(cloud, "alice", "demo", 1)
	serveHeartbeat(cloud, "exp-1")

	entered := make(chan struct{})
	release := make(chan struct{})
	var once sync.Once
	cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		once.Do(func() { close(entered) })
		<-release
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		clienttest.WriteJSON(w, http.StatusCreated, map[string]any{
			"cuid": "exp-1",
			"slug": "slug-drain",
			"name": body["name"],
		})
	})

	mgr := newInstantTestManager()
	s, start, cleanup := newDirectAttachSession(t, mgr, cloud, false, false)
	defer cleanup()

	s.startAttach(attachParams{
		factory:    mgr.factory,
		start:      start,
		resumeMode: api.ResumeNever,
		timing:     fastAttachTiming(),
	})

	<-entered
	// 模拟已经进入收尾排空态后再放行实验创建：接入成功后不应启动后台补传循环。
	if err := s.beginDrain(); err != nil {
		t.Fatalf("beginDrain: %v", err)
	}
	close(release)

	waitAttachDone(t, s)
	if !s.isCloudAttached() {
		t.Fatal("session should attach cloud resources while draining")
	}
	// 验证 s.backfill 从未启动过（stop channel 仍为 nil）。
	if s.backfill == nil {
		t.Fatal("backfill worker should have been created")
	}
}

func TestAttachStopsCleanlyOnRollbackAndShutdown(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, probe.ServerPath, http.StatusOK, map[string]any{"ok": true})
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})

	entered := make(chan struct{})
	release := make(chan struct{})
	var once sync.Once
	cloud.Handle(http.MethodGet, "/api/auth/verify", func(w http.ResponseWriter, r *http.Request) {
		once.Do(func() { close(entered) })
		select {
		case <-r.Context().Done():
			return
		case <-release:
			clienttest.WriteJSON(w, http.StatusOK, map[string]any{"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"})
		}
	})

	mgr := newInstantTestManager()
	s, start, _ := newDirectAttachSession(t, mgr, cloud, false, false)

	s.startAttach(attachParams{
		factory:    mgr.factory,
		start:      start,
		resumeMode: api.ResumeNever,
		timing:     fastAttachTiming(),
	})

	<-entered
	// shutdown -> stopAttach 应取消在途 HTTP 请求并等待接入协程退出。
	ctx, cancel := context.WithTimeout(context.Background(), 2*time.Second)
	defer cancel()
	if err := s.shutdown(ctx); err != nil {
		t.Fatalf("shutdown: %v", err)
	}
	close(release)

	waitManagerCond(t, "attach goroutine exit", func() bool {
		// rollback 同样幂等，确认无残留。
		s.rollback()
		return true
	})

	// 另验证 rollback 场景下也不会泄漏或留下未完成的 run 文件。
	cfg := startSettings(t, cloud)
	cfg.RunDir = filepath.Join(t.TempDir(), "rb")
	cfg.RecoverOnline = true
	if res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
	}); err == nil {
		s2, lookupErr := mgr.session(res.SessionID)
		if lookupErr != nil {
			t.Fatalf("lookup: %v", lookupErr)
		}
		s2.rollback()
	}
}
