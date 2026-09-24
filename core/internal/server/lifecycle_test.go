package server

import (
	"fmt"
	"sync"
	"testing"
	"time"

	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"

	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

// spinup 是测试辅助：把服务置为 READY。
func spinup(t *testing.T, env *testEnv) {
	t.Helper()
	if _, err := env.client.SpinupService(callCtx(t), &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
}

// startRun 是测试辅助：创建一个 run 会话并返回 run_handle。
func startRun(t *testing.T, env *testEnv) string {
	t.Helper()
	resp, err := env.client.DeliverRunStart(
		callCtx(t),
		&corev1.DeliverRunStartRequest{StartRecord: &runv1.StartRecord{}},
	)
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	if !resp.GetSuccess() {
		t.Fatal("DeliverRunStart success = false")
	}
	if resp.GetRunHandle() == "" {
		t.Fatal("DeliverRunStart returned empty run_handle")
	}
	return resp.GetRunHandle()
}

func TestLifecycleStateTransitions(t *testing.T) {
	lc := newLifecycle()
	if lc.Get() != StateNotReady {
		t.Fatalf("initial state = %v, want StateNotReady", lc.Get())
	}
	if !lc.Spinup() {
		t.Fatal("Spinup from NOT_READY must succeed")
	}
	if lc.Get() != StateReady {
		t.Fatalf("state after spinup = %v, want StateReady", lc.Get())
	}
	if !lc.Spinup() {
		t.Fatal("repeated Spinup in READY must succeed (idempotent)")
	}
	lc.BeginStopping()
	if lc.Get() != StateStopping {
		t.Fatalf("state after BeginStopping = %v, want StateStopping", lc.Get())
	}
	lc.BeginStopping() // 幂等
	if lc.Spinup() {
		t.Fatal("Spinup in STOPPING must be rejected")
	}
	lc.Close()
	if lc.Get() != StateClosed {
		t.Fatalf("state after Close = %v, want StateClosed", lc.Get())
	}
	lc.Close() // 幂等
}

func TestSpinupIdempotent(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	spinup(t, env) // 重复 Spinup 在 READY 下幂等成功
}

func TestRunRPCsRejectedBeforeReady(t *testing.T) {
	env := newTestEnv(t)
	ctx := callCtx(t)
	checks := []struct {
		name string
		call func() error
	}{
		{
			"DeliverRunStart",
			func() error {
				_, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{StartRecord: &runv1.StartRecord{}})
				return err
			},
		},
		{
			"UpsertScalars",
			func() error {
				_, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{RunHandle: "any"})
				return err
			},
		},
		{
			"GetOperationStats",
			func() error {
				_, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: "any"})
				return err
			},
		},
		{
			"DeliverRunFinish",
			func() error {
				_, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: "any"})
				return err
			},
		},
		{
			"ConfirmRunFinish",
			func() error {
				_, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: "any"})
				return err
			},
		},
	}
	for _, tc := range checks {
		if err := tc.call(); status.Code(err) != codes.FailedPrecondition {
			t.Fatalf("%s before ready: err = %v, want FailedPrecondition", tc.name, err)
		}
	}
}

func TestUpsertsUnimplementedAfterReady(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env)
	ctx := callCtx(t)
	for _, tc := range []struct {
		name string
		call func() error
	}{
		{
			"UpsertColumns",
			func() error {
				_, err := env.client.UpsertColumns(ctx, &corev1.UpsertColumnsRequest{RunHandle: handle})
				return err
			},
		},
		{
			"UpsertScalars",
			func() error {
				_, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{RunHandle: handle})
				return err
			},
		},
		{
			"UpsertMedia",
			func() error {
				_, err := env.client.UpsertMedia(ctx, &corev1.UpsertMediaRequest{RunHandle: handle})
				return err
			},
		},
		{
			"UpsertLogs",
			func() error {
				_, err := env.client.UpsertLogs(ctx, &corev1.UpsertLogsRequest{RunHandle: handle})
				return err
			},
		},
		{
			"UpsertSaves",
			func() error {
				_, err := env.client.UpsertSaves(ctx, &corev1.UpsertSavesRequest{RunHandle: handle})
				return err
			},
		},
	} {
		if err := tc.call(); status.Code(err) != codes.Unimplemented {
			t.Fatalf("%s after ready: err = %v, want Unimplemented", tc.name, err)
		}
	}
}

func TestDeliverRunStartRequiresStartRecord(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	_, err := env.client.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{})
	if status.Code(err) != codes.InvalidArgument {
		t.Fatalf("DeliverRunStart without start_record: err = %v, want InvalidArgument", err)
	}
}

func TestHandleValidation(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)

	// 空 handle：InvalidArgument
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("empty handle stats err = %v, want InvalidArgument", err)
	}
	// 未知 handle：NotFound
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: "no-such-handle"}); status.Code(err) != codes.NotFound {
		t.Fatalf("unknown handle stats err = %v, want NotFound", err)
	}
}

func TestRunLifecycleCoreStateMapping(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env)
	ctx := callCtx(t)

	// start 后：RUNNING
	resp, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle})
	if err != nil || !resp.GetSuccess() {
		t.Fatalf("stats before finish: err = %v, resp = %v", err, resp)
	}
	if resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("stats state before finish = %v, want RUNNING", resp.GetStats().GetState())
	}

	// finish 后：FINISHED（骨架无数据可排空）
	if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: handle}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	resp, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle})
	if err != nil {
		t.Fatalf("stats after finish: %v", err)
	}
	if resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_FINISHED {
		t.Fatalf("stats state after finish = %v, want FINISHED", resp.GetStats().GetState())
	}

	// confirm：释放会话，不关闭 server
	if _, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); err != nil {
		t.Fatalf("ConfirmRunFinish: %v", err)
	}
	select {
	case <-env.ctrl.Done():
		t.Fatal("ConfirmRunFinish must not shut down the server")
	case <-afterWindow():
	}
	// confirm 后 handle 失效
	if _, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle}); status.Code(err) != codes.NotFound {
		t.Fatalf("stats after confirm err = %v, want NotFound", err)
	}
	if _, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); status.Code(err) != codes.NotFound {
		t.Fatalf("repeat confirm err = %v, want NotFound", err)
	}
}

// TestConfirmBeforeFinishRejected 验证 ConfirmRunFinish 的前置条件：
// finish 未交付时 confirm 返回 FailedPrecondition 且会话保留，
// 之后可完成 finish/confirm 流程。
func TestConfirmBeforeFinishRejected(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env)
	ctx := callCtx(t)

	if _, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); status.Code(err) != codes.FailedPrecondition {
		t.Fatalf("confirm before finish err = %v, want FailedPrecondition", err)
	}
	// 会话未被摘除：stats 报 RUNNING
	resp, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle})
	if err != nil || resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("stats after rejected confirm: err = %v, state = %v, want RUNNING", err, resp.GetStats().GetState())
	}
	// finish 后 confirm 成功，会话被释放
	if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: handle}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if _, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); err != nil {
		t.Fatalf("ConfirmRunFinish after finish: %v", err)
	}
	if _, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle}); status.Code(err) != codes.NotFound {
		t.Fatalf("stats after confirm err = %v, want NotFound", err)
	}
}

// TestFinishConfirmRace 并发执行 finish 与 confirm，验证 registry 的
// 原子校验：confirm 在 finish 标记后才会摘除会话，因此 finish 必然成功；
// confirm 的合法结果为成功（后于 finish）或 FailedPrecondition（先于 finish）。
func TestFinishConfirmRace(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)

	for i := 0; i < 20; i++ {
		handle := startRun(t, env)
		var wg sync.WaitGroup
		var finishErr, confirmErr error
		wg.Add(2)
		go func() {
			defer wg.Done()
			_, finishErr = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: handle})
		}()
		go func() {
			defer wg.Done()
			_, confirmErr = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle})
		}()
		wg.Wait()

		if finishErr != nil {
			t.Fatalf("iter %d: finish must always succeed (session must survive a losing confirm): %v", i, finishErr)
		}
		if confirmErr != nil && status.Code(confirmErr) != codes.FailedPrecondition {
			t.Fatalf("iter %d: confirm err = %v, want nil or FailedPrecondition", i, confirmErr)
		}
		// confirm 输了竞态时会话未释放，补一次 confirm
		if confirmErr != nil {
			if _, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); err != nil {
				t.Fatalf("iter %d: retry confirm after finish: %v", i, err)
			}
		}
		if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle}); status.Code(err) != codes.NotFound {
			t.Fatalf("iter %d: stats after release err = %v, want NotFound", i, err)
		}
	}
}

func TestTwoRunHandleIsolation(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handleA, handleB := startRun(t, env), startRun(t, env)
	if handleA == handleB {
		t.Fatal("two runs must receive distinct run_handles")
	}
	ctx := callCtx(t)

	// 结束 A：A 报 FINISHED，B 保持 RUNNING，互不串扰
	if _, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: handleA}); err != nil {
		t.Fatalf("finish A: %v", err)
	}
	for _, tc := range []struct {
		handle string
		want   operationv1.CoreState
	}{
		{handleA, operationv1.CoreState_CORE_STATE_FINISHED},
		{handleB, operationv1.CoreState_CORE_STATE_RUNNING},
	} {
		resp, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: tc.handle})
		if err != nil {
			t.Fatalf("stats %s: %v", tc.handle, err)
		}
		if resp.GetStats().GetState() != tc.want {
			t.Fatalf("stats %s state = %v, want %v", tc.handle, resp.GetStats().GetState(), tc.want)
		}
	}

	// confirm A：B 可正常走完生命周期
	if _, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handleA}); err != nil {
		t.Fatalf("confirm A: %v", err)
	}
	if _, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: handleB}); err != nil {
		t.Fatalf("finish B: %v", err)
	}
	resp, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handleB})
	if err != nil || resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_FINISHED {
		t.Fatalf("stats B after finish: err = %v, state = %v", err, resp.GetStats().GetState())
	}
}

// TestConcurrentRunSessions 模拟两个并发 client 各自跑完生命周期，
// 配合 -race 验证 session registry 的并发安全与会话隔离。
func TestConcurrentRunSessions(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)

	const clients = 4
	var wg sync.WaitGroup
	errs := make(chan error, clients)
	for i := 0; i < clients; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			ctx := callCtx(t)
			resp, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{StartRecord: &runv1.StartRecord{}})
			if err != nil {
				errs <- fmt.Errorf("client %d start: %w", id, err)
				return
			}
			handle := resp.GetRunHandle()
			for r := 0; r < 20; r++ {
				if _, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle}); err != nil {
					errs <- fmt.Errorf("client %d stats: %w", id, err)
					return
				}
			}
			if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: handle}); err != nil {
				errs <- fmt.Errorf("client %d finish: %w", id, err)
				return
			}
			statsResp, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle})
			if err != nil || statsResp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_FINISHED {
				errs <- fmt.Errorf("client %d stats after finish: err = %v, state = %v", id, err, statsResp.GetStats().GetState())
				return
			}
			if _, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); err != nil {
				errs <- fmt.Errorf("client %d confirm: %w", id, err)
				return
			}
		}(i)
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		t.Error(err)
	}
}

// afterWindow 返回 noShutdownWindow 时长通道，用于断言"未触发关闭"。
func afterWindow() <-chan time.Time {
	return time.After(noShutdownWindow)
}
