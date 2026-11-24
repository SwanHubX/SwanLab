package manager

import (
	"context"
	"encoding/json"
	"errors"
	"net/http"
	"os"
	"sync/atomic"
	"testing"
	"time"

	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager/probe"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func seedLocalRunFile(t *testing.T, cfg config.Settings, name, color string, scalars int, checkpointSeq uint64) {
	t.Helper()
	if err := os.MkdirAll(cfg.RunDir, 0o755); err != nil {
		t.Fatalf("mkdir run_dir: %v", err)
	}
	runFile := runFilePath(t, cfg)
	w, err := store.Create(runFile)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	writeRec := func(rec *recordv1.Record) int64 {
		payload, marshalErr := proto.Marshal(rec)
		if marshalErr != nil {
			t.Fatalf("marshal record: %v", marshalErr)
		}
		if writeErr := w.Write(payload); writeErr != nil {
			t.Fatalf("store.Write: %v", writeErr)
		}
		return w.Offset()
	}
	startRec := &recordv1.Record{
		Num: record.NumStart,
		RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{
			Id:      cfg.RunID,
			Project: "demo",
			Name:    name,
			Color:   color,
		}},
	}
	var seq uint64 = 1
	off := writeRec(startRec)
	var cpOff int64
	if checkpointSeq == seq {
		cpOff = off
	}

	for i := 1; i <= scalars; i++ {
		rec := &recordv1.Record{
			Num: int64(i),
			RecordType: &recordv1.Record_Scalar{Scalar: &datav1.ScalarRecord{
				Key:       "loss",
				Step:      int64(i),
				Type:      columnv1.ColumnType_COLUMN_TYPE_SCALAR,
				Value:     &datav1.ScalarValue{Number: float64(i)},
				Timestamp: timestamppb.New(time.Unix(1_700_000_000, 0)),
			}},
		}
		seq++
		off = writeRec(rec)
		if seq == checkpointSeq {
			cpOff = off
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("store.Close: %v", err)
	}
	if checkpointSeq > 0 {
		if cpOff == 0 {
			t.Fatalf("checkpointSeq %d not reached (wrote %d records)", checkpointSeq, seq)
		}
		if err := store.WriteCheckpoint(runFile+".checkpoint.json", store.Checkpoint{
			Version:  store.CheckpointVersion,
			RunID:    cfg.RunID,
			StoreSeq: checkpointSeq,
			Offset:   cpOff,
		}); err != nil {
			t.Fatalf("WriteCheckpoint: %v", err)
		}
	}
}

func TestDeliverRunStartRecoverOnlineNewExperiment(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var backendUp atomic.Bool
	cloud.Handle(http.MethodGet, probe.ServerPath, func(w http.ResponseWriter, _ *http.Request) {
		if !backendUp.Load() {
			clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]any{"ok": false})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
	serveAuth(cloud)
	serveProject(cloud, "alice", "demo", 1)
	serveHeartbeat(cloud, "exp-recover")
	uploaded := serveMetricsUploadOK(cloud)

	var createBody atomic.Value
	cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		createBody.Store(body)
		clienttest.WriteJSON(w, http.StatusCreated, map[string]any{
			"cuid": "exp-recover",
			"slug": "run-recover-slug",
			"name": body["name"],
		})
	})

	mgr := newInstantTestManager()
	cfg := startSettings(t, cloud)
	cfg.RecoverOnline = true
	cfg.HealthProbeInterval = 0.005
	cfg.HealthProbeTimeout = 0.2
	cfg.HealthDegradeThreshold = 1
	cfg.HealthRecoverThreshold = 1

	startRec := &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "fixed-name",
		Color:   "#112233",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	}
	res, err := mgr.DeliverRunStart(context.Background(), cfg, startRec)
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })

	if res.Path != "" || !res.NewExperiment || res.Name != "fixed-name" {
		t.Fatalf("StartResult = %+v, want empty path, NewExperiment=true, Name=fixed-name", res)
	}
	stats, err := mgr.Stats(res.SessionID)
	if err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if stats.GetCloud().GetState() != operationv1.CloudState_CLOUD_STATE_PENDING {
		t.Fatalf("initial cloud state = %v, want PENDING", stats.GetCloud().GetState())
	}

	// 未接入期间连续写入两条标量，应只写本地而不发往云端。
	for step := int64(1); step <= 2; step++ {
		if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{
			{
				Key:       "loss",
				Step:      step,
				Type:      columnv1.ColumnType_COLUMN_TYPE_SCALAR,
				Value:     &datav1.ScalarValue{Number: float64(step) * 0.1},
				Timestamp: timestamppb.New(time.Unix(1_700_000_000, 0)),
			},
		}); upsertErr != nil {
			t.Fatalf("UpsertScalars: %v", upsertErr)
		}
	}
	if cloud.CountFor(http.MethodPost, "/api/house/metrics") != 0 {
		t.Fatal("records must not be uploaded before cloud attach completes")
	}

	// 后端恢复后自动创建实验并补传全部历史记录。
	backendUp.Store(true)
	s, err := mgr.session(res.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	waitAttachDone(t, s)

	timer := time.NewTimer(5 * time.Second)
	defer timer.Stop()
	select {
	case <-uploaded:
	case <-timer.C:
		t.Fatal("timed out waiting for scalar upload after backend recovery")
	}

	if got := cloud.CountFor(http.MethodPost, "/api/project/alice/demo/experiment"); got != 1 {
		t.Fatalf("experiment create count = %d, want 1", got)
	}
	body, ok := createBody.Load().(map[string]any)
	colors, _ := body["colors"].([]any)
	if !ok || body["name"] != "fixed-name" || body["cuid"] != cfg.RunID ||
		len(colors) != 2 || colors[0] != "#112233" || colors[1] != "#112233" {
		t.Fatalf("create payload = %#v, want fixed cuid/name/colors", createBody.Load())
	}
	statsAfter, err := mgr.Stats(res.SessionID)
	if err != nil {
		t.Fatalf("Stats after attach: %v", err)
	}
	if statsAfter.GetCloud().GetState() != operationv1.CloudState_CLOUD_STATE_ATTACHED ||
		statsAfter.GetCloud().GetPath() != "/alice/demo/run-recover-slug" {
		t.Fatalf("cloud stats after attach = %+v", statsAfter.GetCloud())
	}
}

func TestDeliverRunStartRecoverOnlineWithExistingLocalRun(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var backendUp atomic.Bool
	cloud.Handle(http.MethodGet, probe.ServerPath, func(w http.ResponseWriter, _ *http.Request) {
		if !backendUp.Load() {
			clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]any{"ok": false})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cloud.ServeJSON(http.MethodGet, probe.HousePath, http.StatusOK, map[string]any{"ok": true})
	serveAuth(cloud)
	serveProject(cloud, "alice", "demo", 1)
	serveHeartbeat(cloud, "exp-existing")
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusOK, map[string]any{
		"cuid":      "exp-existing",
		"slug":      "slug-existing",
		"name":      "local-name",
		"createdAt": "2026-01-02T00:00:00.000Z",
	})
	uploaded := serveMetricsUploadOK(cloud)

	mgr := newInstantTestManager()
	cfg := startSettings(t, cloud)
	cfg.RecoverOnline = true
	cfg.HealthProbeInterval = 0.005
	cfg.HealthProbeTimeout = 0.2
	cfg.HealthDegradeThreshold = 1
	cfg.HealthRecoverThreshold = 1

	// 本地已有 Start + 3 条标量，checkpoint 落在第 2 条记录（seq=2），还有 2 条待补传。
	seedLocalRunFile(t, cfg, "local-name", "#654321", 3, 2)

	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Resume:  runv1.ResumeMode_RESUME_MODE_ALLOW,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })

	if res.NewExperiment {
		t.Fatal("NewExperiment = true, want false when reusing existing local run")
	}
	if res.Run.GetName() != "local-name" || res.Run.GetColor() != "#654321" {
		t.Fatalf("StartRecord lost local identity: name=%q color=%q", res.Run.GetName(), res.Run.GetColor())
	}
	stats, err := mgr.Stats(res.SessionID)
	if err != nil {
		t.Fatalf("Stats: %v", err)
	}
	if stats.GetTotalRecords() != 2 {
		t.Fatalf("TotalRecords = %d, want 2 unconfirmed records from checkpoint", stats.GetTotalRecords())
	}

	backendUp.Store(true)
	s, err := mgr.session(res.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	waitAttachDone(t, s)

	timer := time.NewTimer(5 * time.Second)
	defer timer.Stop()
	select {
	case <-uploaded:
	case <-timer.C:
		t.Fatal("timed out waiting for checkpointed replay after attach")
	}

	// 恢复模式接入不拉取云端摘要接口。
	if summaryCalls := cloud.CountFor(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-existing"); summaryCalls != 0 {
		t.Fatalf("ExperimentSummary called %d times, want 0", summaryCalls)
	}
	// 确认只补传 checkpoint 之后的记录（step 2 与 step 3，不重复传 step 1）。
	waitManagerCond(t, "uploaded count reaches 2", func() bool {
		return s.tracker.Snapshot().GetUploadedRecords() == 2
	})
}

func TestDeliverRunStartRecoverOnlineRejectsResumeWithoutLocalFile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	mgr := newInstantTestManager()

	for _, mode := range []runv1.ResumeMode{runv1.ResumeMode_RESUME_MODE_ALLOW, runv1.ResumeMode_RESUME_MODE_MUST} {
		cfg := startSettings(t, cloud)
		cfg.RecoverOnline = true
		_, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
			Id:      cfg.RunID,
			Project: "demo",
			Resume:  mode,
		})
		if !errors.Is(err, ErrInvalidStartRecord) {
			t.Fatalf("mode %v got err=%v, want ErrInvalidStartRecord", mode, err)
		}
	}
}

func TestDeliverRunStartRecoverOnlineStoreDegradedDisablesSuspension(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	mgr := newInstantTestManager()
	s, _, cleanup := newDirectAttachSession(t, mgr, cloud, true, false)
	defer cleanup()

	if !s.shouldSuspendUploads() {
		t.Fatal("unattached online session with healthy local store should suspend uploads")
	}
	s.storeDegraded.Store(true)
	if s.shouldSuspendUploads() {
		t.Fatal("store degraded must disable upload suspension even before cloud attach")
	}
}
