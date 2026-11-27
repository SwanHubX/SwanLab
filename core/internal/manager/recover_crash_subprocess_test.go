// 本文件用真实子进程验证进程被强制终止后的恢复：
// 进程没有执行 Shutdown，也不会刷盘或写 checkpoint，恢复只能依赖已经保存的 run 文件和 checkpoint。
package manager

import (
	"bytes"
	"context"
	"encoding/json"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/store"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

const (
	crashHelperEnv   = "SWANLAB_TEST_CRASH_HELPER"
	crashStageEnv    = "SWANLAB_TEST_CRASH_STAGE"
	crashAPIHostEnv  = "SWANLAB_TEST_CRASH_API_HOST"
	crashRunIDEnv    = "SWANLAB_TEST_CRASH_RUN_ID"
	crashRunDirEnv   = "SWANLAB_TEST_CRASH_RUN_DIR"
	crashResumeEnv   = "SWANLAB_TEST_CRASH_RESUME"
	crashBeforeStage = "before_experiment_created"
	crashAfterStage  = "after_experiment_created"
	crashMidStage    = "mid_backfill_checkpoint"
)

// TestRecoverCrashHelperProcess 只作为父测试的子进程运行；普通测试直接跳过。
func TestRecoverCrashHelperProcess(t *testing.T) {
	if os.Getenv(crashHelperEnv) != "1" {
		t.Skip("subprocess crash helper")
	}
	cfg := config.Defaults()
	cfg.Mode = config.ModeOnline
	cfg.RunID = os.Getenv(crashRunIDEnv)
	cfg.RunDir = os.Getenv(crashRunDirEnv)
	cfg.APIHost = os.Getenv(crashAPIHostEnv)
	cfg.WebHost = cfg.APIHost
	cfg.APIKey = "key-alice"
	cfg.RecordBatch = 1
	cfg.RecordInterval = 0.01
	cfg.RecoverOnline = true
	cfg.HealthProbeInterval = 0.005
	cfg.HealthProbeTimeout = 0.2
	cfg.HealthDegradeThreshold = 1
	cfg.HealthRecoverThreshold = 1

	resume := runv1.ResumeMode_RESUME_MODE_NEVER
	if raw := os.Getenv(crashResumeEnv); raw != "" {
		if mode, ok := runv1.ResumeMode_value[raw]; ok {
			resume = runv1.ResumeMode(mode)
		}
	}

	mgr := newInstantTestManager()
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "crash-" + os.Getenv(crashStageEnv),
		Resume:  resume,
	})
	if err != nil {
		t.Fatalf("helper DeliverRunStart: %v", err)
	}

	if os.Getenv(crashStageEnv) != crashMidStage {
		if err := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{
			scalarAtStep("loss", 1, 0.5),
		}); err != nil {
			t.Fatalf("helper UpsertScalars: %v", err)
		}
	}

	// 父进程会在目标阶段强制结束本进程；这里保持运行，避免执行正常的 Shutdown。
	select {}
}

type lockedBuffer struct {
	mu sync.Mutex
	b  bytes.Buffer
}

func (l *lockedBuffer) Write(p []byte) (int, error) {
	l.mu.Lock()
	defer l.mu.Unlock()
	return l.b.Write(p)
}

func (l *lockedBuffer) String() string {
	l.mu.Lock()
	defer l.mu.Unlock()
	return l.b.String()
}

type crashChild struct {
	cmd  *exec.Cmd
	out  *lockedBuffer
	done chan error
}

func startCrashChild(t *testing.T, cfg config.Settings, stage string, resume runv1.ResumeMode) *crashChild {
	t.Helper()
	cmd := exec.Command(os.Args[0], "-test.run=^TestRecoverCrashHelperProcess$")
	out := &lockedBuffer{}
	cmd.Stdout = out
	cmd.Stderr = out
	cmd.Env = append(os.Environ(),
		crashHelperEnv+"=1",
		crashStageEnv+"="+stage,
		crashAPIHostEnv+"="+cfg.APIHost,
		crashRunIDEnv+"="+cfg.RunID,
		crashRunDirEnv+"="+cfg.RunDir,
		crashResumeEnv+"="+resume.String(),
	)
	if err := cmd.Start(); err != nil {
		t.Fatalf("start crash helper: %v", err)
	}
	child := &crashChild{cmd: cmd, out: out, done: make(chan error, 1)}
	go func() { child.done <- cmd.Wait() }()
	t.Cleanup(func() {
		if cmd.ProcessState == nil {
			_ = cmd.Process.Kill()
			select {
			case err := <-child.done:
				_ = err
			case <-time.After(5 * time.Second):
			}
		}
	})
	return child
}

func waitCrashCondition(t *testing.T, child *crashChild, desc string, cond func() bool) {
	t.Helper()
	deadline := time.Now().Add(15 * time.Second)
	for !cond() {
		if child != nil {
			select {
			case err := <-child.done:
				t.Fatalf("crash helper exited before %s: %v\n%s", desc, err, child.out.String())
			default:
			}
		}
		if time.Now().After(deadline) {
			output := ""
			if child != nil {
				output = child.out.String()
			}
			t.Fatalf("timed out waiting for %s\n%s", desc, output)
		}
		time.Sleep(10 * time.Millisecond)
	}
}

func killCrashChild(t *testing.T, child *crashChild) {
	t.Helper()
	if err := child.cmd.Process.Kill(); err != nil {
		t.Fatalf("kill crash helper: %v", err)
	}
	select {
	case err := <-child.done:
		_ = err
	case <-time.After(5 * time.Second):
		t.Fatal("crash helper did not exit after kill")
	}
}

func waitRunFileFlushed(t *testing.T, cfg config.Settings) {
	t.Helper()
	waitCrashCondition(t, nil, "run file flush", func() bool {
		matches, _ := filepath.Glob(filepath.Join(cfg.RunDir, "run-*.swanlab"))
		if len(matches) != 1 {
			return false
		}
		info, err := os.Stat(matches[0])
		return err == nil && info.Size() > store.HeaderLen
	})
}

func waitCheckpointSeq(t *testing.T, cfg config.Settings, want int64) {
	t.Helper()
	path := runFilePath(t, cfg) + ".checkpoint.json"
	deadline := time.Now().Add(15 * time.Second)
	for {
		data, err := os.ReadFile(path)
		if err == nil {
			var cp store.Checkpoint
			if json.Unmarshal(data, &cp) == nil && cp.StoreSeq >= want {
				return
			}
		}
		if time.Now().After(deadline) {
			t.Fatalf("timed out waiting for checkpoint %s store_seq >= %d: %v", path, want, err)
		}
		time.Sleep(20 * time.Millisecond)
	}
}

type crashCloud struct {
	cloud     *clienttest.Cloud
	backendUp atomic.Bool

	mu          sync.Mutex
	createdAt   map[string]int
	experiments map[string]string
	steps       map[string][]int64
	stopCalls   map[string]int
	stopStates  map[string]string

	blocked chan string
	release chan struct{}
}

func newCrashCloud(t *testing.T, runID string, block func(runID string, steps []int64) bool) *crashCloud {
	t.Helper()
	cc := &crashCloud{
		cloud:       clienttest.NewCloud(t),
		createdAt:   make(map[string]int),
		experiments: make(map[string]string),
		steps:       make(map[string][]int64),
		stopCalls:   make(map[string]int),
		stopStates:  make(map[string]string),
		blocked:     make(chan string, 1),
		release:     make(chan struct{}),
	}
	t.Cleanup(func() { close(cc.release) })

	cc.cloud.ServeHealth(&cc.backendUp, nil, http.StatusBadGateway)
	cc.cloud.Handle(http.MethodGet, "/api/auth/verify", func(w http.ResponseWriter, _ *http.Request) {
		if !cc.backendUp.Load() {
			clienttest.WriteJSON(w, http.StatusBadGateway, map[string]any{"message": "backend down"})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"uid": 1, "username": "alice", "createdAt": "2026-01-01T00:00:00Z"})
	})
	cc.cloud.Handle(http.MethodGet, "/api/project/alice/demo", func(w http.ResponseWriter, _ *http.Request) {
		if !cc.backendUp.Load() {
			clienttest.WriteJSON(w, http.StatusBadGateway, map[string]any{"message": "backend down"})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{
			"cuid": "proj-1", "name": "demo", "path": "/@alice/demo", "version": 1,
			"_count": map[string]any{"experiments": 1},
			"group":  map[string]any{"username": "alice"},
		})
	})
	cc.cloud.Handle(http.MethodPost, "/api/project/alice/demo/experiment", func(w http.ResponseWriter, r *http.Request) {
		if !cc.backendUp.Load() {
			clienttest.WriteJSON(w, http.StatusBadGateway, map[string]any{"message": "backend down"})
			return
		}
		var body struct {
			Cuid string `json:"cuid"`
			Name string `json:"name"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		cuid := body.Cuid
		if cuid == "" {
			cuid = "exp-generated"
		}
		cc.mu.Lock()
		status := http.StatusOK
		if _, exists := cc.experiments[cuid]; !exists {
			cc.experiments[cuid] = cuid
			cc.createdAt[cuid]++
			status = http.StatusCreated
		}
		cc.mu.Unlock()
		clienttest.WriteJSON(w, status, map[string]any{
			"cuid": cuid, "slug": cuid, "name": body.Name, "createdAt": "2026-01-02T00:00:00.000Z",
		})
	})
	cc.cloud.ServeSummary("", runID, nil)
	cc.cloud.ServeHeartbeat(runID)
	cc.cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, r *http.Request) {
		var payload struct {
			ExperimentID string `json:"experimentId"`
			Type         string `json:"type"`
			Metrics      []struct {
				Index int64 `json:"index"`
			} `json:"metrics"`
		}
		_ = json.NewDecoder(r.Body).Decode(&payload)
		if payload.Type != "scalar" {
			clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
			return
		}
		indices := make([]int64, 0, len(payload.Metrics))
		for _, item := range payload.Metrics {
			indices = append(indices, item.Index)
		}
		if block != nil && block(payload.ExperimentID, indices) {
			select {
			case cc.blocked <- payload.ExperimentID:
			default:
			}
			select {
			case <-cc.release:
			case <-r.Context().Done():
			}
			return
		}
		cc.mu.Lock()
		cc.steps[payload.ExperimentID] = append(cc.steps[payload.ExperimentID], indices...)
		cc.mu.Unlock()
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cc.cloud.Handle(http.MethodPut, "/api/project/alice/demo/runs/"+runID+"/state", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			State string `json:"state"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		cc.mu.Lock()
		cc.stopCalls[runID]++
		cc.stopStates[runID] = body.State
		cc.mu.Unlock()
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cc.cloud.ServeSeries("alice", "demo")
	return cc
}

func (c *crashCloud) uploadedSteps(runID string) []int64 {
	c.mu.Lock()
	defer c.mu.Unlock()
	out := append([]int64(nil), c.steps[runID]...)
	sort.Slice(out, func(i, j int) bool { return out[i] < out[j] })
	return out
}

func (c *crashCloud) experimentCount(runID string) int {
	c.mu.Lock()
	defer c.mu.Unlock()
	return c.createdAt[runID]
}

func (c *crashCloud) stopState(runID string) (string, int) {
	c.mu.Lock()
	defer c.mu.Unlock()
	return c.stopStates[runID], c.stopCalls[runID]
}

func crashSettings(t *testing.T, cc *crashCloud, runID string) config.Settings {
	t.Helper()
	cfg := config.Defaults()
	cfg.Mode = config.ModeOnline
	cfg.RunID = runID
	cfg.RunDir = filepath.Join(t.TempDir(), "run")
	cfg.RecordBatch = 1
	cfg.RecordInterval = 0.01
	cfg.APIHost = cc.cloud.URL()
	cfg.WebHost = cc.cloud.URL()
	cfg.APIKey = "key-alice"
	cfg.RecoverOnline = true
	cfg.HealthProbeInterval = 0.005
	cfg.HealthProbeTimeout = 0.2
	cfg.HealthDegradeThreshold = 1
	cfg.HealthRecoverThreshold = 1
	return cfg
}

func resumeAfterCrash(t *testing.T, cfg config.Settings, addStep int64) {
	t.Helper()
	cfg.RecoverOnline = false
	mgr := newInstantTestManager()
	res, err := mgr.DeliverRunStart(context.Background(), cfg, &runv1.StartRecord{
		Id:      cfg.RunID,
		Project: "demo",
		Name:    "restarted",
		Resume:  runv1.ResumeMode_RESUME_MODE_ALLOW,
	})
	if err != nil {
		t.Fatalf("restart DeliverRunStart: %v", err)
	}
	t.Cleanup(func() { closeStartedSession(t, mgr, res.SessionID) })
	if addStep > 0 {
		if upsertErr := mgr.UpsertScalars(context.Background(), res.SessionID, []*datav1.ScalarRecord{
			scalarAtStep("loss", addStep, 0.25),
		}); upsertErr != nil {
			t.Fatalf("restart UpsertScalars: %v", upsertErr)
		}
	}
	if finishErr := mgr.DeliverRunFinish(context.Background(), res.SessionID, finishRecord(runv1.RunState_RUN_STATE_FINISHED, "")); finishErr != nil {
		t.Fatalf("restart DeliverRunFinish: %v", finishErr)
	}
	conf, err := mgr.Confirm(context.Background(), res.SessionID)
	if err != nil || conf.ReleaseFailed {
		t.Fatalf("restart Confirm = (%+v, %v)", conf, err)
	}
}

func assertCrashResult(t *testing.T, cc *crashCloud, runID string, wantSteps []int64) {
	t.Helper()
	if got := cc.experimentCount(runID); got != 1 {
		t.Fatalf("created experiments for %s = %d, want 1", runID, got)
	}
	got := cc.uploadedSteps(runID)
	if len(got) != len(wantSteps) {
		t.Fatalf("uploaded steps = %v, want %v", got, wantSteps)
	}
	for i := range wantSteps {
		if got[i] != wantSteps[i] {
			t.Fatalf("uploaded steps = %v, want %v", got, wantSteps)
		}
	}
	state, calls := cc.stopState(runID)
	if calls != 1 || state != "FINISHED" {
		t.Fatalf("stop state = (%q, %d), want (FINISHED, 1)", state, calls)
	}
}

func TestRecoverOnlineProcessKillStages(t *testing.T) {
	t.Run(crashBeforeStage, func(t *testing.T) {
		runID := "kill-before-exp"
		cc := newCrashCloud(t, runID, nil)
		cfg := crashSettings(t, cc, runID)
		child := startCrashChild(t, cfg, crashBeforeStage, runv1.ResumeMode_RESUME_MODE_NEVER)
		waitRunFileFlushed(t, cfg)
		killCrashChild(t, child)
		if got := cc.experimentCount(runID); got != 0 {
			t.Fatalf("experiment created before backend recovery = %d, want 0", got)
		}

		cc.backendUp.Store(true)
		resumeAfterCrash(t, cfg, 2)
		assertCrashResult(t, cc, runID, []int64{1, 2})
	})

	t.Run(crashAfterStage, func(t *testing.T) {
		runID := "kill-after-exp"
		var blockedOnce atomic.Bool
		cc := newCrashCloud(t, runID, func(expID string, _ []int64) bool {
			if expID != runID {
				return false
			}
			return blockedOnce.CompareAndSwap(false, true)
		})
		cc.backendUp.Store(true)
		cfg := crashSettings(t, cc, runID)
		child := startCrashChild(t, cfg, crashAfterStage, runv1.ResumeMode_RESUME_MODE_NEVER)
		waitCrashCondition(t, child, "experiment created", func() bool {
			return cc.experimentCount(runID) == 1
		})
		select {
		case <-cc.blocked:
		case <-time.After(15 * time.Second):
			t.Fatalf("timed out waiting for first metrics upload\n%s", child.out.String())
		}
		killCrashChild(t, child)

		resumeAfterCrash(t, cfg, 2)
		assertCrashResult(t, cc, runID, []int64{1, 2})
	})

	t.Run(crashMidStage, func(t *testing.T) {
		runID := "kill-mid-backfill"
		var blockBackfill atomic.Bool
		blockBackfill.Store(true)
		cc := newCrashCloud(t, runID, func(expID string, steps []int64) bool {
			if expID != runID || !blockBackfill.Load() {
				return false
			}
			for _, step := range steps {
				if step >= 5 {
					return true
				}
			}
			return false
		})
		cc.backendUp.Store(true)
		cfg := crashSettings(t, cc, runID)
		seedLocalRunFile(t, cfg, "local-mid", "#112233", 6, 4)
		child := startCrashChild(t, cfg, crashMidStage, runv1.ResumeMode_RESUME_MODE_ALLOW)
		select {
		case <-cc.blocked:
		case <-time.After(15 * time.Second):
			t.Fatalf("timed out waiting for blocked backfill upload\n%s", child.out.String())
		}
		// 等待 step 4 上传成功并把 checkpoint 推进到 seq 5，然后强制结束进程；
		// 重启后应从 checkpoint 之后继续，而不是从第一条记录重传。
		waitCheckpointSeq(t, cfg, 5)
		killCrashChild(t, child)
		blockBackfill.Store(false)

		resumeAfterCrash(t, cfg, 0)
		assertCrashResult(t, cc, runID, []int64{4, 5, 6})
	})
}
