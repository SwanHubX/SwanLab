package manager

import (
	"context"
	"encoding/json"
	"errors"
	"net/http"
	"os"
	"path/filepath"
	"reflect"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/cron/health"
	"github.com/swanhubx/swanlab/core/internal/manager/probe"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func serveAuth(cloud *clienttest.Cloud) {
	cloud.ServeJSON(http.MethodGet, "/api/auth/verify", http.StatusOK, map[string]any{
		"uid":       1,
		"username":  "alice",
		"createdAt": "2026-01-01T00:00:00Z",
	})
}

func serveHeartbeat(cloud *clienttest.Cloud, experimentID string) {
	cloud.ServeJSON(http.MethodPost, "/api/house/experiments/"+experimentID+"/heartbeat",
		http.StatusOK, map[string]any{"ok": true})
}

func serveProject(cloud *clienttest.Cloud, owner, name string, experiments int) {
	cloud.ServeJSON(http.MethodGet, "/api/project/"+owner+"/"+name, http.StatusOK, map[string]any{
		"cuid":       "proj-1",
		"name":       name,
		"path":       "/" + owner + "/" + name,
		"visibility": "PRIVATE",
		"version":    2,
		"_count":     map[string]any{"experiments": experiments},
		"group":      map[string]any{"username": owner},
	})
}

func serveNewExperiment(cloud *clienttest.Cloud, owner, project, cuid, slug string) {
	cloud.Handle(http.MethodPost, "/api/project/"+owner+"/"+project+"/experiment", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			Name string `json:"name"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		clienttest.WriteJSON(w, http.StatusCreated, map[string]any{
			"cuid": cuid,
			"slug": slug,
			"name": body.Name,
		})
	})
}

func startSettings(t *testing.T, cloud *clienttest.Cloud) config.Settings {
	t.Helper()
	settings := config.Defaults()
	settings.RunID = "run-1"
	settings.RunDir = filepath.Join(t.TempDir(), "run")
	settings.RecordBatch = 16
	settings.RecordInterval = 0.01
	settings.APIKey = "smoke-key"
	settings.APIHost = cloud.URL()
	settings.WebHost = cloud.URL()
	return settings
}

func runFilePath(t *testing.T, settings config.Settings) string {
	t.Helper()
	run, err := runctx.New(settings)
	if err != nil {
		t.Fatalf("runctx.New: %v", err)
	}
	path, err := run.RunFile()
	if err != nil {
		t.Fatalf("RunFile: %v", err)
	}
	return path
}

func closeStartedSession(t *testing.T, mgr *Manager, sessionID string) {
	t.Helper()
	s, err := mgr.session(sessionID)
	if err != nil {
		return
	}
	if s.cancel != nil {
		s.cancel()
	}
	if s.heartbeat != nil {
		s.heartbeat.Stop(context.Background())
	}
	if s.watcher != nil {
		s.watcher.Stop()
	}
	if s.uploader != nil {
		_ = s.uploader.Close()
	}
	if s.intake != nil {
		ctx, cancel := context.WithTimeout(context.Background(), time.Second)
		_ = s.intake.Close(ctx)
		cancel()
	}
	if s.store != nil {
		_ = s.store.Close()
	}
	s.closeHealth()
}

func TestDeliverRunStartLocalWritesStore(t *testing.T) {
	settings := config.Defaults()
	settings.Mode = config.ModeLocal
	settings.RunID = "run-local"
	settings.RunDir = filepath.Join(t.TempDir(), "run")
	settings.RecordBatch = 16
	settings.RecordInterval = 0.01

	mgr := newTestManager()
	record := &runv1.StartRecord{Project: "demo", Id: "run-local"}
	result, err := mgr.DeliverRunStart(context.Background(), settings, record)
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	defer closeStartedSession(t, mgr, result.SessionID)
	if result.SessionID == "" {
		t.Fatal("local start must return a session id")
	}
	if result.Path != "" || result.GlobalStep != -1 || result.GlobalSystemStep != -1 {
		t.Fatalf("local start result = %+v, want empty path and step -1", result)
	}

	scalar := &datav1.ScalarRecord{
		Key:   "loss",
		Step:  1,
		Type:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		Value: &datav1.ScalarValue{Number: 0.5},
	}
	if upsertErr := mgr.UpsertScalars(context.Background(), result.SessionID, []*datav1.ScalarRecord{scalar}); upsertErr != nil {
		t.Fatalf("UpsertScalars: %v", upsertErr)
	}
	if finishErr := mgr.DeliverRunFinish(context.Background(), result.SessionID, &runv1.FinishRecord{
		State: runv1.RunState_RUN_STATE_FINISHED,
	}); finishErr != nil {
		t.Fatalf("DeliverRunFinish: %v", finishErr)
	}
	confirm, err := mgr.Confirm(context.Background(), result.SessionID)
	if err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if confirm.ReleaseFailed {
		t.Fatalf("Confirm release failed: %s", confirm.Message)
	}
	if _, sessionErr := mgr.session(result.SessionID); sessionErr == nil {
		t.Fatal("local session must be released after confirm")
	}
	info, err := os.Stat(runFilePath(t, settings))
	if err != nil {
		t.Fatalf("stat run file: %v", err)
	}
	if info.Size() < store.HeaderLen {
		t.Fatalf("run file size = %d, want >= %d", info.Size(), store.HeaderLen)
	}
}

func TestDeliverRunStartNewExperiment(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 7)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-slug")

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	record := &runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"}

	result, err := mgr.DeliverRunStart(context.Background(), settings, record)
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, result.SessionID) })

	if want := builder.GenerateName(7); result.Name != want {
		t.Fatalf("name = %q, want %q", result.Name, want)
	}
	if want := builder.GenerateColor(7); result.Run.Color != want {
		t.Fatalf("color = %q, want %q", result.Run.Color, want)
	}
	if result.Run.Name != result.Name {
		t.Fatalf("run name = %q, want %q", result.Run.Name, result.Name)
	}
	if result.Path != "/alice/demo/run-slug" {
		t.Fatalf("path = %q, want /alice/demo/run-slug", result.Path)
	}
	if !result.NewExperiment || result.GlobalStep != -1 || result.GlobalSystemStep != -1 {
		t.Fatalf("unexpected start result: %+v", result)
	}
	if result.Run.Workspace != "alice" || result.Run.Project != "demo" {
		t.Fatalf("start record not filled: %+v", result.Run)
	}
	mustState(t, mgr, result.SessionID, RunStateRunning)

	runFile := runFilePath(t, settings)
	info, err := os.Stat(runFile)
	if err != nil {
		t.Fatalf("run file missing: %v", err)
	}
	if info.Size() < store.HeaderLen {
		t.Fatalf("run file size = %d, want >= %d", info.Size(), store.HeaderLen)
	}

	requests := cloud.RequestsFor(http.MethodPost, "/api/project/alice/demo/experiment")
	if len(requests) != 1 {
		t.Fatalf("create requests = %d, want 1", len(requests))
	}
	var body map[string]any
	if err := json.Unmarshal(requests[0].Body, &body); err != nil {
		t.Fatalf("decode create body: %v", err)
	}
	if body["cuid"] != "run-1" {
		t.Fatalf("create body cuid = %v, want run-1", body["cuid"])
	}
	if body["name"] != result.Name {
		t.Fatalf("create body name = %v, want %q", body["name"], result.Name)
	}
}

func TestDeliverRunStartResume(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 3)
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusOK, map[string]any{
		"cuid":      "exp-1",
		"slug":      "run-1",
		"name":      "resume",
		"createdAt": "2026-01-02T00:00:00.000Z",
	})
	cloud.ServeJSON(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-1", http.StatusOK, map[string]any{
		"scalar": []map[string]any{{"key": "loss", "step": 12}},
		"media":  []map[string]any{{"key": "image", "step": 3}},
		"log":    []map[string]any{{"key": "console", "step": 21}},
	})

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	record := &runv1.StartRecord{
		Project:   "demo",
		Workspace: "alice",
		Id:        "run-1",
		Name:      "resume",
		Color:     "#123456",
		Resume:    runv1.ResumeMode_RESUME_MODE_ALLOW,
	}

	result, err := mgr.DeliverRunStart(context.Background(), settings, record)
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, result.SessionID) })

	if result.NewExperiment {
		t.Fatal("200 response must be reported as an existing experiment")
	}
	if result.GlobalStep != 12 {
		t.Fatalf("global step = %d, want 12", result.GlobalStep)
	}
	if result.GlobalSystemStep != -1 {
		t.Fatalf("global system step = %d, want -1", result.GlobalSystemStep)
	}
	if result.Name != "resume" {
		t.Fatalf("response name = %q, want experiment name", result.Name)
	}
	if result.Run.Color != "#123456" {
		t.Fatalf("explicit color was overwritten: %q", result.Run.Color)
	}
	if cloud.CountFor(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-1") != 1 {
		t.Fatal("resume must fetch the experiment summary once")
	}
	s, err := mgr.session(result.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	if !s.metrics.Defined("loss") {
		t.Fatal("resume metrics were not attached to the session")
	}
	if minStep, ok := s.metrics.MinStep("loss"); !ok || minStep != 12 {
		t.Fatalf("loss min step = %d (defined=%v), want 12", minStep, ok)
	}
	if s.epoch == nil || s.epoch.Value() != 21 {
		t.Fatalf("console epoch not restored: %+v", s.epoch)
	}
	if s.counter == nil || s.counter.Value() != 0 {
		t.Fatalf("record counter must start at 0: %+v", s.counter)
	}
}

func TestDeliverRunStartUsesGroupWorkspace(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	cloud.ServeJSON(http.MethodGet, "/api/project/alice/demo", http.StatusOK, map[string]any{
		"cuid":    "proj-1",
		"name":    "demo",
		"path":    "/org/demo",
		"version": 1,
		"_count":  map[string]any{"experiments": 1},
		"group":   map[string]any{"username": "org"},
	})
	serveNewExperiment(cloud, "org", "demo", "exp-1", "run-1")

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	record := &runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"}

	result, err := mgr.DeliverRunStart(context.Background(), settings, record)
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, result.SessionID) })

	if result.Run.Workspace != "org" {
		t.Fatalf("workspace = %q, want org", result.Run.Workspace)
	}
	if result.Path != "/org/demo/run-1" {
		t.Fatalf("path = %q, want /org/demo/run-1", result.Path)
	}
	if cloud.CountFor(http.MethodPost, "/api/project/org/demo/experiment") != 1 {
		t.Fatal("experiment must be created under the project group owner")
	}
}

func TestDeliverRunStartSkipStore(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 1)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-1")

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	settings.SkipStore = true
	record := &runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"}

	result, err := mgr.DeliverRunStart(context.Background(), settings, record)
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, result.SessionID) })

	if _, err = os.Stat(settings.RunDir); !os.IsNotExist(err) {
		t.Fatalf("skip_store must not create run dir: err=%v", err)
	}
	s, err := mgr.session(result.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	if s.store != nil {
		t.Fatal("skip_store must not create the local store")
	}
}

func TestDeliverRunStartRollbackRemovesRunFile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 1)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-1")

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	settings.SaveSize = 0 // runctx 不校验，NewSave 在 store 创建后失败
	record := &runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"}

	if _, err := mgr.DeliverRunStart(context.Background(), settings, record); err == nil {
		t.Fatal("DeliverRunStart must fail when save limits are invalid")
	}
	if _, err := os.Stat(runFilePath(t, settings)); !os.IsNotExist(err) {
		t.Fatalf("rollback must remove the run file: err=%v", err)
	}
	mgr.mu.RLock()
	count := len(mgr.sessions)
	mgr.mu.RUnlock()
	if count != 0 {
		t.Fatalf("sessions = %d, want 0", count)
	}
}

func TestDeliverRunStartKeepsExistingRunFile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 1)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-1")

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	if err := os.MkdirAll(settings.RunDir, 0o755); err != nil {
		t.Fatal(err)
	}
	runFile := runFilePath(t, settings)
	if err := os.WriteFile(runFile, []byte("keep"), 0o644); err != nil {
		t.Fatal(err)
	}
	record := &runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"}

	if _, err := mgr.DeliverRunStart(context.Background(), settings, record); err == nil {
		t.Fatal("DeliverRunStart must fail when the run file already exists")
	}
	data, err := os.ReadFile(runFile)
	if err != nil {
		t.Fatalf("existing run file was removed: %v", err)
	}
	if string(data) != "keep" {
		t.Fatalf("existing run file changed: %q", data)
	}
}

func TestDeliverRunStartRejectsInvalidInput(t *testing.T) {
	mgr := newTestManager()
	settings := config.Defaults()
	settings.RunID = "run-1"
	settings.RunDir = t.TempDir()

	if _, err := mgr.DeliverRunStart(context.Background(), settings, nil); err == nil {
		t.Fatal("nil start record must be rejected")
	}
	record := &runv1.StartRecord{Project: "demo"}
	record.Resume = runv1.ResumeMode(99)
	if _, err := mgr.DeliverRunStart(context.Background(), settings, record); err == nil {
		t.Fatal("unknown resume mode must be rejected")
	}
	record.Resume = runv1.ResumeMode_RESUME_MODE_NEVER
	settings.RunDir = ""
	if _, err := mgr.DeliverRunStart(context.Background(), settings, record); err == nil {
		t.Fatal("empty run_dir must be rejected")
	}
}

func TestSummaryFromAPI(t *testing.T) {
	got := summaryFromAPI(&api.Summary{
		Log:    []api.ColumnStep{{Key: "console", Step: 3}},
		Media:  []api.ColumnStep{{Key: "image", Step: 7}},
		Scalar: []api.ColumnStep{{Key: "loss", Step: 9}, {Key: "__swanlab__.gpu", Step: 11}},
	})
	want := &metrics.Summary{
		Log:    []metrics.SummaryItem{{Key: "console", Step: 3}},
		Media:  []metrics.SummaryItem{{Key: "image", Step: 7}},
		Scalar: []metrics.SummaryItem{{Key: "loss", Step: 9}, {Key: "__swanlab__.gpu", Step: 11}},
	}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("summary = %#v, want %#v", got, want)
	}
}

func TestSummaryFromAPINilAndEmpty(t *testing.T) {
	if got := summaryFromAPI(nil); got != nil {
		t.Fatalf("nil summary = %#v, want nil", got)
	}
	if got := summaryFromAPI(&api.Summary{}); !reflect.DeepEqual(got, &metrics.Summary{}) {
		t.Fatalf("empty summary = %#v, want empty", got)
	}
}

func TestResumeModeFromProto(t *testing.T) {
	cases := []struct {
		in   runv1.ResumeMode
		want api.ResumeMode
	}{
		{runv1.ResumeMode_RESUME_MODE_NEVER, api.ResumeNever},
		{runv1.ResumeMode_RESUME_MODE_ALLOW, api.ResumeAllow},
		{runv1.ResumeMode_RESUME_MODE_MUST, api.ResumeMust},
	}
	for _, tc := range cases {
		got, err := resumeModeFromProto(tc.in)
		if err != nil {
			t.Fatalf("resumeModeFromProto(%v): %v", tc.in, err)
		}
		if got != tc.want {
			t.Fatalf("resumeModeFromProto(%v) = %q, want %q", tc.in, got, tc.want)
		}
	}
	if _, err := resumeModeFromProto(runv1.ResumeMode(99)); err == nil {
		t.Fatal("unknown resume mode must return an error")
	}
}

func TestDeliverRunStartFailsFastWhenControlUnreachable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	settings := startSettings(t, cloud)
	cloud.Close()

	mgr := newTestManager()
	_, err := mgr.DeliverRunStart(context.Background(), settings,
		&runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"})
	if !errors.Is(err, ErrServerUnreachable) {
		t.Fatalf("err = %v, want ErrServerUnreachable", err)
	}
	mgr.mu.RLock()
	sessions := len(mgr.sessions)
	mgr.mu.RUnlock()
	if sessions != 0 {
		t.Fatalf("sessions = %d, want 0", sessions)
	}
	if got := mgr.healthReg.Len(); got != 0 {
		t.Fatalf("sites = %d, want 0", got)
	}
}

func TestDeliverRunStartNormalizesHostBeforeProbe(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 7)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-slug")

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	settings.APIHost = cloud.URL() + "/ignored/path?ignored=true"
	result, err := mgr.DeliverRunStart(context.Background(), settings,
		&runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"})
	if err != nil {
		t.Fatalf("DeliverRunStart with path and query: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, result.SessionID) })
	if got := mgr.healthReg.Len(); got != 1 {
		t.Fatalf("sites = %d, want 1", got)
	}
}

func TestDeliverRunStartDegradedWhenDataPlaneDown(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveAuth(cloud)
	serveHeartbeat(cloud, "exp-1")
	serveProject(cloud, "alice", "demo", 7)
	serveNewExperiment(cloud, "alice", "demo", "exp-1", "run-slug")
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusServiceUnavailable, nil)

	mgr := newTestManager()
	settings := startSettings(t, cloud)
	result, err := mgr.DeliverRunStart(context.Background(), settings,
		&runv1.StartRecord{Project: "demo", Workspace: "alice", Id: "run-1"})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, result.SessionID) })

	mustState(t, mgr, result.SessionID, RunStateRunning)
	s, err := mgr.session(result.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	snap := s.site.Snapshot()
	if snap.State != health.StateDegraded || snap.Reason != health.ReasonProbeHouse {
		t.Fatalf("snapshot = %+v, want degraded by swanlab-house", snap)
	}
	if !s.uploader.Suspended() {
		t.Fatal("uploader should start suspended when swanlab-house is down")
	}
	if got := mgr.healthReg.Len(); got != 1 {
		t.Fatalf("sites = %d, want 1", got)
	}
}
