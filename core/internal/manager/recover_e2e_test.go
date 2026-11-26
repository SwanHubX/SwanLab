package manager

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

func scalarAtStep(key string, step int64, val float64) *datav1.ScalarRecord {
	return &datav1.ScalarRecord{
		Key:       key,
		Step:      step,
		Type:      columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		Value:     &datav1.ScalarValue{Number: val},
		Timestamp: timestamppb.New(time.Unix(1_700_000_000+step, 0)),
	}
}

func fastRecoverSettings(t *testing.T, cloud *clienttest.Cloud, runID, runDir string) config.Settings {
	t.Helper()
	cfg := startSettings(t, cloud)
	cfg.RunID = runID
	cfg.RunDir = runDir
	cfg.RecoverOnline = true
	cfg.HealthProbeInterval = 0.005
	cfg.HealthProbeTimeout = 0.2
	cfg.HealthDegradeThreshold = 1
	cfg.HealthRecoverThreshold = 1
	return cfg
}

func waitScalarSteps(t *testing.T, mu *sync.Mutex, got *[]int64, wantCount int) []int64 {
	t.Helper()
	waitManagerCond(t, fmt.Sprintf("%d scalar steps uploaded", wantCount), func() bool {
		mu.Lock()
		defer mu.Unlock()
		return len(*got) >= wantCount
	})
	mu.Lock()
	defer mu.Unlock()
	out := append([]int64(nil), (*got)...)
	sort.Slice(out, func(i, j int) bool { return out[i] < out[j] })
	return out
}

func TestRecoverOnlineMidAttachOutageAndPostAttachOutage(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var serverUp, houseUp atomic.Bool

	cloud.ServeRun(clienttest.RunSpec{
		ExpID:      "exp-outage",
		Slug:       "slug-outage",
		ServerUp:   &serverUp,
		HouseUp:    &houseUp,
		DownStatus: http.StatusBadGateway,
	})
	cloud.Handle(http.MethodGet, "/api/auth/verify", func(w http.ResponseWriter, _ *http.Request) {
		if !serverUp.Load() {
			clienttest.WriteJSON(w, http.StatusBadGateway, map[string]any{"message": "bad gateway"})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"})
	})

	var expCreates atomic.Int32
	cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		if !serverUp.Load() {
			clienttest.WriteJSON(w, http.StatusBadGateway, map[string]any{"message": "bad gateway"})
			return
		}
		expCreates.Add(1)
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		clienttest.WriteJSON(w, http.StatusCreated, map[string]any{
			"cuid": "exp-outage",
			"slug": "slug-outage",
			"name": body["name"],
		})
	})

	var mu sync.Mutex
	var uploadedSteps []int64
	var stopStates []string
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, r *http.Request) {
		if !houseUp.Load() {
			clienttest.WriteJSON(w, http.StatusBadGateway, map[string]any{"message": "house down"})
			return
		}
		var payload struct {
			Type    string `json:"type"`
			Metrics []struct {
				Index int64 `json:"index"`
			} `json:"metrics"`
		}
		_ = json.NewDecoder(r.Body).Decode(&payload)
		if payload.Type == "scalar" {
			mu.Lock()
			for _, item := range payload.Metrics {
				uploadedSteps = append(uploadedSteps, item.Index)
			}
			mu.Unlock()
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cloud.Handle(http.MethodPut, "/api/project/alice/demo/runs/exp-outage/state", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			State string `json:"state"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		mu.Lock()
		stopStates = append(stopStates, body.State)
		mu.Unlock()
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	mgr := newInstantTestManager()
	cfg := fastRecoverSettings(t, cloud, "run-outage", filepath.Join(t.TempDir(), "run"))
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "outage-flow",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })

	s, err := mgr.session(res.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}

	// 1. 云端实验创建前写入第 1 条指标，此时服务端不可用。
	if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 1, 1.0)}); upsertErr != nil {
		t.Fatalf("UpsertScalars step 1: %v", upsertErr)
	}

	// 2. 主站先恢复但指标服务仍返回 502：只创建 1 次实验，上传保持暂停等待。
	serverUp.Store(true)
	waitAttachDone(t, s)
	if got := expCreates.Load(); got != 1 {
		t.Fatalf("expCreates after attach = %d, want 1", got)
	}
	waitManagerCond(t, "uploader suspended while house down", func() bool {
		up := s.uploader()
		return up != nil && up.Suspended()
	})
	mu.Lock()
	preHouseRecoverCount := len(uploadedSteps)
	mu.Unlock()
	if preHouseRecoverCount != 0 {
		t.Fatalf("uploadedSteps while house down = %d, want 0", preHouseRecoverCount)
	}

	// 3. 指标服务恢复后上传第 1 条，并验证创建实验后写入的第 2 条实时上传。
	houseUp.Store(true)
	waitScalarSteps(t, &mu, &uploadedSteps, 1)
	if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 2, 0.5)}); upsertErr != nil {
		t.Fatalf("UpsertScalars step 2: %v", upsertErr)
	}
	waitScalarSteps(t, &mu, &uploadedSteps, 2)

	// 4. 创建实验后再次遇到指标服务宕机，写入第 3、4 条；恢复后结束运行，应全部上传完成且不重建实验。
	houseUp.Store(false)
	waitManagerCond(t, "house probe degrades", func() bool {
		up := s.uploader()
		return up != nil && up.Suspended()
	})
	if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{
		scalarAtStep("loss", 3, 0.25),
		scalarAtStep("loss", 4, 0.125),
	}); upsertErr != nil {
		t.Fatalf("UpsertScalars steps 3..4: %v", upsertErr)
	}

	houseUp.Store(true)
	if finErr := mgr.DeliverRunFinish(context.Background(), res.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
		t.Fatalf("DeliverRunFinish: %v", finErr)
	}
	conf, err := mgr.Confirm(context.Background(), res.SessionID)
	if err != nil {
		t.Fatalf("Confirm: %v", err)
	}
	if conf.ReleaseFailed {
		t.Fatalf("Confirm ReleaseFailed=true, msg=%q", conf.Message)
	}
	if got := expCreates.Load(); got != 1 {
		t.Fatalf("final expCreates = %d, want 1", got)
	}
	gotSteps := waitScalarSteps(t, &mu, &uploadedSteps, 4)
	for i := int64(1); i <= 4; i++ {
		if gotSteps[i-1] != i {
			t.Fatalf("uploadedSteps = %v, want [1 2 3 4]", gotSteps)
		}
	}
	mu.Lock()
	defer mu.Unlock()
	if len(stopStates) != 1 || stopStates[0] != "FINISHED" {
		t.Fatalf("stopStates = %v, want [FINISHED]", stopStates)
	}
}

func TestRecoverOnlineConfirmCanceledThenRetried(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var backendUp atomic.Bool
	cloud.ServeRun(clienttest.RunSpec{
		ExpID:      "exp-cancel",
		Slug:       "slug-cancel",
		ServerUp:   &backendUp,
		DownStatus: http.StatusBadGateway,
	})
	uploaded := serveMetricsUploadOK(cloud)

	mgr := newInstantTestManager()
	cfg := fastRecoverSettings(t, cloud, "run-cancel", filepath.Join(t.TempDir(), "run"))
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "cancel-retry",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })

	if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 1, 0.9)}); upsertErr != nil {
		t.Fatalf("UpsertScalars: %v", upsertErr)
	}
	if finErr := mgr.DeliverRunFinish(context.Background(), res.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
		t.Fatalf("DeliverRunFinish: %v", finErr)
	}

	// 第一次 Confirm 超时取消后，会话仍保留以待重试。
	shortCtx, cancel := context.WithTimeout(context.Background(), 40*time.Millisecond)
	defer cancel()
	if _, confErr := mgr.Confirm(shortCtx, res.SessionID); !errors.Is(confErr, context.DeadlineExceeded) {
		t.Fatalf("canceled Confirm err = %v, want context.DeadlineExceeded", confErr)
	}
	if mgr.SessionCount() != 1 {
		t.Fatalf("SessionCount after canceled confirm = %d, want 1", mgr.SessionCount())
	}

	backendUp.Store(true)
	conf, err := mgr.Confirm(context.Background(), res.SessionID)
	if err != nil {
		t.Fatalf("retried Confirm: %v", err)
	}
	if conf.ReleaseFailed || mgr.SessionCount() != 0 {
		t.Fatalf("retried Confirm = %+v, SessionCount = %d", conf, mgr.SessionCount())
	}
	select {
	case <-uploaded:
	case <-time.After(2 * time.Second):
		t.Fatal("expected buffered scalar to upload on retried Confirm")
	}
}

func TestRecoverOnlineZeroRecordFinishWaitsUntilAttach(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var backendUp atomic.Bool
	cloud.ServeRun(clienttest.RunSpec{
		ExpID:      "exp-zero",
		Slug:       "slug-zero",
		ServerUp:   &backendUp,
		DownStatus: http.StatusBadGateway,
	})

	mgr := newInstantTestManager()
	cfg := fastRecoverSettings(t, cloud, "run-zero", filepath.Join(t.TempDir(), "run"))
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "zero-records",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })

	if finErr := mgr.DeliverRunFinish(context.Background(), res.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
		t.Fatalf("DeliverRunFinish: %v", finErr)
	}

	done := make(chan error, 1)
	go func() {
		conf, confErr := mgr.Confirm(context.Background(), res.SessionID)
		if confErr != nil {
			done <- confErr
			return
		}
		if conf.ReleaseFailed {
			done <- errors.New(conf.Message)
			return
		}
		done <- nil
	}()

	// 即使没有指标记录，也不能因为上传队列为空而提前放行 Confirm。
	select {
	case earlyErr := <-done:
		t.Fatalf("Confirm returned before cloud attach: %v", earlyErr)
	case <-time.After(60 * time.Millisecond):
	}

	backendUp.Store(true)
	select {
	case waitErr := <-done:
		if waitErr != nil {
			t.Fatalf("Confirm after attach failed: %v", waitErr)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("Confirm timed out after backend recovery")
	}
	if got := cloud.CountFor(http.MethodPost, "/api/project/alice/demo/experiment"); got != 1 {
		t.Fatalf("experiment creates = %d, want 1", got)
	}
}

func TestRecoverOnlineHighVolumeBacklogMemoryAndIntegrity(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var down atomic.Bool
	cloud.ServeHealth(&down, &down, http.StatusBadGateway)

	mgr := newInstantTestManager()
	cfg := fastRecoverSettings(t, cloud, "run-stress", filepath.Join(t.TempDir(), "run"))
	cfg.RecordBatch = 1000
	cfg.RecordInterval = 0.005

	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "stress-200k",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}

	runtime.GC()
	var memBefore, memAfter runtime.MemStats
	runtime.ReadMemStats(&memBefore)

	const (
		batchSize    = 1000
		totalBatches = 200
		totalScalars = batchSize * totalBatches
	)
	batch := make([]*datav1.ScalarRecord, batchSize)
	for b := range totalBatches {
		for i := range batchSize {
			step := int64(b*batchSize + i + 1)
			batch[i] = scalarAtStep("loss", step, 1.0/float64(step))
		}
		if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, batch); upsertErr != nil {
			t.Fatalf("UpsertScalars batch %d: %v", b, upsertErr)
		}
	}

	if finErr := mgr.DeliverRunFinish(context.Background(), res.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
		t.Fatalf("DeliverRunFinish: %v", finErr)
	}

	runtime.GC()
	runtime.ReadMemStats(&memAfter)
	if memAfter.HeapAlloc > memBefore.HeapAlloc {
		delta := memAfter.HeapAlloc - memBefore.HeapAlloc
		const maxAllowedHeapDelta = 256 << 20
		if delta > maxAllowedHeapDelta {
			t.Fatalf("heap grew by %d bytes (> %d), expected bounded memory while waiting", delta, maxAllowedHeapDelta)
		}
	}

	shutdownCtx, cancel := context.WithTimeout(context.Background(), 10*time.Second)
	defer cancel()
	if shutErr := mgr.Shutdown(shutdownCtx); shutErr != nil {
		t.Fatalf("Shutdown: %v", shutErr)
	}

	runFile := runFilePath(t, cfg)
	r, err := store.Open(runFile)
	if err != nil {
		t.Fatalf("store.Open: %v", err)
	}
	defer func() { _ = r.Close() }()

	var scalarCount int
	var hasStart, hasFinish bool
	for {
		payload, scanErr := r.Scan()
		if errors.Is(scanErr, io.EOF) {
			break
		}
		if scanErr != nil {
			t.Fatalf("r.Scan: %v", scanErr)
		}
		if payload == nil {
			break
		}
		var rec recordv1.Record
		if unmarshalErr := proto.Unmarshal(payload, &rec); unmarshalErr != nil {
			t.Fatalf("proto.Unmarshal: %v", unmarshalErr)
		}
		switch rec.GetRecordType().(type) {
		case *recordv1.Record_Start:
			hasStart = true
		case *recordv1.Record_Scalar:
			scalarCount++
		case *recordv1.Record_Finish:
			hasFinish = true
		}
	}
	if !hasStart || !hasFinish || scalarCount != totalScalars {
		t.Fatalf("persisted records: start=%v finish=%v scalars=%d, want true/true/%d", hasStart, hasFinish, scalarCount, totalScalars)
	}
}

func setupIdempotentRecoveryCloud(t *testing.T, backendUp *atomic.Bool) (*clienttest.Cloud, *[]int64, *sync.Mutex, map[string]int) {
	t.Helper()
	cloud := clienttest.NewCloud(t)
	cloud.ServeHealth(backendUp, nil, http.StatusBadGateway)
	cloud.ServeAuth("")
	cloud.ServeProject(clienttest.ProjectSpec{})
	cloud.ServeSeries("alice", "demo")

	var mu sync.Mutex
	var steps []int64
	createdByCuid := make(map[string]int)

	cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			Cuid string `json:"cuid"`
			Name string `json:"name"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		cuid := body.Cuid
		if cuid == "" {
			cuid = "exp-generated"
		}
		mu.Lock()
		createdByCuid[cuid]++
		count := createdByCuid[cuid]
		mu.Unlock()
		status := http.StatusCreated
		if count > 1 {
			status = http.StatusOK
		}
		clienttest.WriteJSON(w, status, map[string]any{
			"cuid":      cuid,
			"slug":      cuid,
			"name":      body.Name,
			"createdAt": "2026-01-02T00:00:00.000Z",
		})
	})
	cloud.Handle(http.MethodGet, "/api/project/alice/demo/runs/crash-after", func(w http.ResponseWriter, _ *http.Request) {
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{
			"cuid":      "crash-after",
			"slug":      "crash-after",
			"name":      "crash-after-name",
			"createdAt": "2026-01-02T00:00:00.000Z",
		})
	})
	cloud.ServeSummary("", "crash-after", nil)
	cloud.ServeHeartbeat("crash-pre")
	cloud.ServeHeartbeat("crash-after")
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, r *http.Request) {
		var payload struct {
			Type    string `json:"type"`
			Metrics []struct {
				Index int64 `json:"index"`
			} `json:"metrics"`
		}
		_ = json.NewDecoder(r.Body).Decode(&payload)
		if payload.Type == "scalar" {
			mu.Lock()
			for _, item := range payload.Metrics {
				steps = append(steps, item.Index)
			}
			mu.Unlock()
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cloud.ServeRunState("alice", "demo", "crash-pre")
	cloud.ServeRunState("alice", "demo", "crash-after")
	return cloud, &steps, &mu, createdByCuid
}

func runShutdownBeforeExperimentCreated(t *testing.T) {
	var backendUp atomic.Bool
	cloud, steps, mu, createdByCuid := setupIdempotentRecoveryCloud(t, &backendUp)
	runDir := filepath.Join(t.TempDir(), "crash-pre-dir")

	mgr1 := newInstantTestManager()
	cfg1 := fastRecoverSettings(t, cloud, "crash-pre", runDir)
	res1, err := mgr1.DeliverRunStart(context.Background(), cfg1, &runv1.StartRecord{
		Id:      cfg1.RunID,
		Project: "demo",
		Name:    "crash-pre-name",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("mgr1 DeliverRunStart: %v", err)
	}
	if upsertErr := mgr1.UpsertScalars(context.Background(), res1.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 1, 1.0)}); upsertErr != nil {
		t.Fatalf("mgr1 UpsertScalars: %v", upsertErr)
	}
	waitManagerCond(t, "scalar 1 flushed before crash", func() bool {
		matches, _ := filepath.Glob(filepath.Join(runDir, "run-*.swanlab"))
		if len(matches) != 1 {
			return false
		}
		info, statErr := os.Stat(matches[0])
		return statErr == nil && info.Size() > 0
	})
	if shutErr := mgr1.Shutdown(context.Background()); shutErr != nil {
		t.Fatalf("mgr1 Shutdown: %v", shutErr)
	}

	// Shutdown 后重新启动：允许续传并恢复网络，应上传中断前的记录，并且只创建一个实验。
	backendUp.Store(true)
	mgr2 := newInstantTestManager()
	cfg2 := fastRecoverSettings(t, cloud, "crash-pre", runDir)
	res2, err := mgr2.DeliverRunStart(context.Background(), cfg2, &runv1.StartRecord{
		Id:      cfg2.RunID,
		Project: "demo",
		Resume:  runv1.ResumeMode_RESUME_MODE_ALLOW,
	})
	if err != nil {
		t.Fatalf("mgr2 DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr2, res2.SessionID) })
	if upsertErr := mgr2.UpsertScalars(context.Background(), res2.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 2, 0.5)}); upsertErr != nil {
		t.Fatalf("mgr2 UpsertScalars: %v", upsertErr)
	}
	if finErr := mgr2.DeliverRunFinish(context.Background(), res2.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
		t.Fatalf("mgr2 DeliverRunFinish: %v", finErr)
	}
	conf, err := mgr2.Confirm(context.Background(), res2.SessionID)
	if err != nil || conf.ReleaseFailed {
		t.Fatalf("mgr2 Confirm = (%+v, %v)", conf, err)
	}
	got := waitScalarSteps(t, mu, steps, 2)
	if len(got) != 2 || got[0] != 1 || got[1] != 2 {
		t.Fatalf("steps = %v, want [1 2]", got)
	}
	mu.Lock()
	defer mu.Unlock()
	if len(createdByCuid) != 1 || createdByCuid["crash-pre"] != 1 {
		t.Fatalf("createdByCuid = %v, want 1 unique experiment", createdByCuid)
	}
}

func runShutdownAfterExperimentCreated(t *testing.T) {
	var backendUp atomic.Bool
	backendUp.Store(true)
	cloud, _, mu, createdByCuid := setupIdempotentRecoveryCloud(t, &backendUp)
	runDir := filepath.Join(t.TempDir(), "crash-after-dir")

	mgr1 := newInstantTestManager()
	cfg1 := fastRecoverSettings(t, cloud, "crash-after", runDir)
	res1, err := mgr1.DeliverRunStart(context.Background(), cfg1, &runv1.StartRecord{
		Id:      cfg1.RunID,
		Project: "demo",
		Name:    "crash-after-name",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("mgr1 DeliverRunStart: %v", err)
	}
	s1, err := mgr1.session(res1.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	waitAttachDone(t, s1)
	if upsertErr := mgr1.UpsertScalars(context.Background(), res1.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 1, 1.0)}); upsertErr != nil {
		t.Fatalf("UpsertScalars: %v", upsertErr)
	}
	if shutErr := mgr1.Shutdown(context.Background()); shutErr != nil {
		t.Fatalf("mgr1 Shutdown: %v", shutErr)
	}

	// 以同一 cuid 重拉恢复时复用同一实验，不产生第二个实验。
	mgr2 := newInstantTestManager()
	cfg2 := fastRecoverSettings(t, cloud, "crash-after", runDir)
	res2, err := mgr2.DeliverRunStart(context.Background(), cfg2, &runv1.StartRecord{
		Id:      cfg2.RunID,
		Project: "demo",
		Resume:  runv1.ResumeMode_RESUME_MODE_ALLOW,
	})
	if err != nil {
		t.Fatalf("mgr2 DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr2, res2.SessionID) })
	if finErr := mgr2.DeliverRunFinish(context.Background(), res2.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
		t.Fatalf("mgr2 DeliverRunFinish: %v", finErr)
	}
	conf, err := mgr2.Confirm(context.Background(), res2.SessionID)
	if err != nil || conf.ReleaseFailed {
		t.Fatalf("mgr2 Confirm = (%+v, %v)", conf, err)
	}
	mu.Lock()
	defer mu.Unlock()
	if len(createdByCuid) != 1 || createdByCuid["crash-after"] != 2 {
		t.Fatalf("createdByCuid = %v, want single experiment id reused", createdByCuid)
	}
}

func TestRecoverOnlineShutdownStagesNoDuplicateExperiment(t *testing.T) {
	t.Run("shutdown_before_experiment_created", runShutdownBeforeExperimentCreated)
	t.Run("shutdown_after_experiment_created", runShutdownAfterExperimentCreated)
}

// TestRecoverOnlineFinishFailureDuringAttachKeepsCloudResources 覆盖“本地结束失败”
// 和“云端实验正在创建”同时发生的场景：实验已经创建时，连接任务仍要保存云端资源；
// Confirm 随后完成上传并上报终态，不能因为本地结束失败就误报云端永久失败。
func TestRecoverOnlineFinishFailureDuringAttachKeepsCloudResources(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeRun(clienttest.RunSpec{ExpID: "exp-finish-race", Slug: "slug-finish-race"})
	uploaded := serveMetricsUploadOK(cloud)

	experimentStarted := make(chan struct{})
	releaseExperiment := make(chan struct{})
	var releaseOnce sync.Once
	release := func() { releaseOnce.Do(func() { close(releaseExperiment) }) }
	t.Cleanup(release)
	var startedOnce sync.Once
	var stopCalls atomic.Int32
	var stopState atomic.Value
	cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		startedOnce.Do(func() { close(experimentStarted) })
		<-releaseExperiment
		var body map[string]any
		_ = json.NewDecoder(r.Body).Decode(&body)
		clienttest.WriteJSON(w, http.StatusCreated, map[string]any{
			"cuid": "exp-finish-race",
			"slug": "slug-finish-race",
			"name": body["name"],
		})
	})
	cloud.Handle(http.MethodPut, "/api/project/alice/demo/runs/exp-finish-race/state", func(w http.ResponseWriter, r *http.Request) {
		stopCalls.Add(1)
		var body struct {
			State string `json:"state"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		stopState.Store(body.State)
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	mgr := newInstantTestManager()
	cfg := fastRecoverSettings(t, cloud, "run-finish-race", filepath.Join(t.TempDir(), "run"))
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "finish-race",
		Resume:  runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })

	s, err := mgr.session(res.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}

	// 1. 阻塞实验创建请求，让连接任务停在“实验已请求、资源尚未保存”的位置。
	select {
	case <-experimentStarted:
	case <-time.After(5 * time.Second):
		t.Fatal("timed out waiting for experiment create request")
	}
	if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{
		scalarAtStep("loss", 1, 0.5),
	}); upsertErr != nil {
		t.Fatalf("UpsertScalars: %v", upsertErr)
	}

	// 2. 模拟 DeliverRunFinish 写入 FinishRecord 失败：会话进入 FAILED，FinishRecord 尚未写入。
	//    这里直接构造相同状态，避免依赖队列是否有空位来触发失败；Confirm 会继续完成后续步骤。
	if started, drainErr := s.beginDrainOwned(); drainErr != nil || !started {
		t.Fatalf("beginDrainOwned = (%v, %v), want (true, nil)", started, drainErr)
	}
	s.backfill.Stop()
	s.watcher.Stop()
	s.prepareFinish(finishRecord(runv1.RunState_RUN_STATE_FINISHED, ""))
	finishFailure := errors.New("simulated local finish commit failure")
	s.fail(finishFailure)
	originalFailure := s.failure()
	if !errors.Is(originalFailure, finishFailure) {
		t.Fatalf("finish failure = %v, want wrapped local commit failure", originalFailure)
	}

	// 3. 放开实验创建请求：FAILED 会话仍要保存云端资源，而不是关闭上传器并误报永久失败。
	release()
	waitAttachDone(t, s)
	if state, path, attachErr := s.cloudStatus(); state != operationv1.CloudState_CLOUD_STATE_ATTACHED ||
		path != "/alice/demo/slug-finish-race" || attachErr != nil {
		t.Fatalf("cloudStatus = (%v, %q, %v), want ATTACHED", state, path, attachErr)
	}
	if failure := s.failure(); !errors.Is(failure, finishFailure) || strings.Contains(failure.Error(), "cloud attach") {
		t.Fatalf("finish failure was overwritten: %v", failure)
	}

	// 4. 第二次 Confirm 完成本地结束、上传积压记录、等待上传完成并上报终态。
	conf, err := mgr.Confirm(context.Background(), res.SessionID)
	if err != nil || conf.ReleaseFailed {
		t.Fatalf("Confirm = (%+v, %v)", conf, err)
	}
	timer := time.NewTimer(5 * time.Second)
	defer timer.Stop()
	select {
	case <-uploaded:
	case <-timer.C:
		t.Fatal("timed out waiting for backfill after resumed finish")
	}
	if got := stopCalls.Load(); got != 1 {
		t.Fatalf("stop_experiment calls = %d, want 1", got)
	}
	if got, _ := stopState.Load().(string); got != "FINISHED" {
		t.Fatalf("stop state = %q, want FINISHED", got)
	}
}

func TestRecoverOnlineMultiSessionSharedSiteIsolation(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var backendUp atomic.Bool
	cloud.ServeHealth(&backendUp, nil, http.StatusBadGateway)
	cloud.Handle(http.MethodGet, "/api/auth/verify", func(w http.ResponseWriter, r *http.Request) {
		switch r.Header.Get("Authorization") {
		case "ApiKey key-alice":
			clienttest.WriteJSON(w, http.StatusOK, map[string]any{"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"})
		case "ApiKey key-bob":
			clienttest.WriteJSON(w, http.StatusOK, map[string]any{"uid": 2, "username": "bob", "createdAt": "2026-01-01T00:00:00Z"})
		default:
			clienttest.WriteJSON(w, http.StatusUnauthorized, map[string]any{"message": "unknown key"})
		}
	})
	for _, spec := range []struct {
		owner, proj, expID, slug string
	}{
		{"alice", "proj-a", "exp-alice", "slug-alice"},
		{"bob", "proj-b", "exp-bob", "slug-bob"},
	} {
		cloud.ServeProject(clienttest.ProjectSpec{Owner: spec.owner, Name: spec.proj})
		cloud.ServeExperiment(clienttest.ExperimentSpec{Owner: spec.owner, Project: spec.proj, Cuid: spec.expID, Slug: spec.slug})
		cloud.ServeHeartbeat(spec.expID)
		cloud.ServeSeries(spec.owner, spec.proj)
		cloud.ServeRunState(spec.owner, spec.proj, spec.expID)
	}

	var mu sync.Mutex
	expSteps := map[string][]int64{}
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, r *http.Request) {
		var payload struct {
			ExperimentID string `json:"experimentId"`
			Type         string `json:"type"`
			Metrics      []struct {
				Index int64 `json:"index"`
			} `json:"metrics"`
		}
		_ = json.NewDecoder(r.Body).Decode(&payload)
		if payload.Type == "scalar" {
			mu.Lock()
			for _, item := range payload.Metrics {
				expSteps[payload.ExperimentID] = append(expSteps[payload.ExperimentID], item.Index)
			}
			mu.Unlock()
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	mgr := newInstantTestManager()
	cfgA := fastRecoverSettings(t, cloud, "iso-a", filepath.Join(t.TempDir(), "a"))
	cfgA.APIKey = "key-alice"
	cfgB := fastRecoverSettings(t, cloud, "iso-b", filepath.Join(t.TempDir(), "b"))
	cfgB.APIKey = "key-bob"

	resA, err := mgr.DeliverRunStart(context.Background(), cfgA, &runv1.StartRecord{
		Id: cfgA.RunID, Project: "proj-a", Name: "run-a", Resume: runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart A: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, resA.SessionID) })

	resB, err := mgr.DeliverRunStart(context.Background(), cfgB, &runv1.StartRecord{
		Id: cfgB.RunID, Project: "proj-b", Name: "run-b", Resume: runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart B: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, resB.SessionID) })

	if upsertErr := mgr.UpsertScalars(context.Background(), resA.SessionID, []*datav1.ScalarRecord{scalarAtStep("acc", 11, 0.9)}); upsertErr != nil {
		t.Fatalf("UpsertScalars A: %v", upsertErr)
	}
	if upsertErr := mgr.UpsertScalars(context.Background(), resB.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 22, 0.1)}); upsertErr != nil {
		t.Fatalf("UpsertScalars B: %v", upsertErr)
	}

	backendUp.Store(true)
	for _, sid := range []string{resA.SessionID, resB.SessionID} {
		if finErr := mgr.DeliverRunFinish(context.Background(), sid, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
			t.Fatalf("DeliverRunFinish %s: %v", sid, finErr)
		}
		conf, confErr := mgr.Confirm(context.Background(), sid)
		if confErr != nil || conf.ReleaseFailed {
			t.Fatalf("Confirm %s = (%+v, %v)", sid, conf, confErr)
		}
	}

	if mgr.SessionCount() != 0 || mgr.healthReg.Len() != 0 {
		t.Fatalf("leaked resources: sessions=%d healthSites=%d", mgr.SessionCount(), mgr.healthReg.Len())
	}
	mu.Lock()
	defer mu.Unlock()
	if got := expSteps["exp-alice"]; len(got) != 1 || got[0] != 11 {
		t.Fatalf("exp-alice steps = %v, want [11]", got)
	}
	if got := expSteps["exp-bob"]; len(got) != 1 || got[0] != 22 {
		t.Fatalf("exp-bob steps = %v, want [22]", got)
	}
}

func TestRecoverOnlinePermanentAuthFailureCleansProbeAndSession(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeHealth(nil, nil, 0)
	var authCalls atomic.Int32
	cloud.Handle(http.MethodGet, "/api/auth/verify", func(w http.ResponseWriter, _ *http.Request) {
		if authCalls.Add(1) == 1 {
			clienttest.WriteJSON(w, http.StatusBadGateway, map[string]any{"message": "temporary outage"})
			return
		}
		clienttest.WriteJSON(w, http.StatusUnauthorized, map[string]any{"code": "auth_error", "message": "expired key"})
	})

	mgr := newInstantTestManager()
	cfg := fastRecoverSettings(t, cloud, "run-auth-perm", filepath.Join(t.TempDir(), "run"))
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id: cfg.RunID, Project: "demo", Name: "auth-perm", Resume: runv1.ResumeMode_RESUME_MODE_NEVER,
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })

	s, err := mgr.session(res.SessionID)
	if err != nil {
		t.Fatalf("session: %v", err)
	}
	waitAttachDone(t, s)
	if state, _, _ := s.cloudStatus(); state != operationv1.CloudState_CLOUD_STATE_FAILED {
		t.Fatalf("cloud state = %v, want FAILED", state)
	}

	// 永久失败后训练写入仍应继续落盘成功。
	if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{scalarAtStep("loss", 1, 0.5)}); upsertErr != nil {
		t.Fatalf("UpsertScalars after permanent failure: %v", upsertErr)
	}
	if finErr := mgr.DeliverRunFinish(context.Background(), res.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finErr != nil {
		t.Fatalf("DeliverRunFinish: %v", finErr)
	}

	first, err := mgr.Confirm(context.Background(), res.SessionID)
	if err != nil || !first.ReleaseFailed || first.Retryable || !first.Released || !strings.Contains(first.Message, "expired key") {
		t.Fatalf("first Confirm = (%+v, %v)", first, err)
	}
	if mgr.SessionCount() != 0 || mgr.healthReg.Len() != 0 {
		t.Fatalf("after cleanup: sessions=%d healthSites=%d, want 0/0", mgr.SessionCount(), mgr.healthReg.Len())
	}
	second, err := mgr.Confirm(context.Background(), res.SessionID)
	if err != nil || second != first {
		t.Fatalf("second Confirm = (%+v, %v), want identical %+v", second, err, first)
	}
}
