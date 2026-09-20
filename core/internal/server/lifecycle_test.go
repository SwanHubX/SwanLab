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

// spinup 是测试辅助：用正确的 owner/auth token 把服务置为 READY。
func spinup(t *testing.T, env *testEnv) {
	t.Helper()
	if _, err := env.client.SpinupService(authCtx(t, testAuthToken), &corev1.SpinupServiceRequest{OwnerToken: testOwnerToken}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
}

// startRun 是测试辅助：创建一个 run 会话并返回 run_handle。
func startRun(t *testing.T, env *testEnv) string {
	t.Helper()
	resp, err := env.client.DeliverRunStart(
		authCtx(t, testAuthToken),
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

func TestAuthInterceptor(t *testing.T) {
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	ctx := plainCtx(t)

	// 无 auth metadata：被 interceptor 拒绝
	if _, err := env.client.SpinupService(ctx, &corev1.SpinupServiceRequest{OwnerToken: testOwnerToken}); status.Code(err) != codes.PermissionDenied {
		t.Fatalf("no-auth SpinupService err = %v, want PermissionDenied", err)
	}
	// 错误 auth token
	wrong := authCtx(t, "wrong-auth")
	if _, err := env.client.SpinupService(wrong, &corev1.SpinupServiceRequest{OwnerToken: testOwnerToken}); status.Code(err) != codes.PermissionDenied {
		t.Fatalf("wrong-auth SpinupService err = %v, want PermissionDenied", err)
	}
	select {
	case <-env.ctrl.Done():
		t.Fatal("auth failure must not trigger shutdown")
	case <-afterWindow():
	}
}

func TestSpinupOwnerCheck(t *testing.T) {
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	// 正确 auth 但错误 owner token
	_, err := env.client.SpinupService(authCtx(t, testAuthToken), &corev1.SpinupServiceRequest{OwnerToken: "wrong-owner"})
	if status.Code(err) != codes.PermissionDenied {
		t.Fatalf("wrong-owner SpinupService err = %v, want PermissionDenied", err)
	}
	// 未配置 owner token 的服务一律拒绝
	env2 := newTestEnv(t, "", testAuthToken)
	_, err = env2.client.SpinupService(authCtx(t, testAuthToken), &corev1.SpinupServiceRequest{OwnerToken: "any"})
	if status.Code(err) != codes.PermissionDenied {
		t.Fatalf("empty-configured-owner SpinupService err = %v, want PermissionDenied", err)
	}
}

func TestSpinupIdempotent(t *testing.T) {
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	spinup(t, env)
	spinup(t, env) // 重复 Spinup 在 READY 下幂等成功
}

func TestRunRPCsRejectedBeforeReady(t *testing.T) {
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	ctx := authCtx(t, testAuthToken)
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
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	spinup(t, env)
	handle := startRun(t, env)
	ctx := authCtx(t, testAuthToken)
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
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	spinup(t, env)
	_, err := env.client.DeliverRunStart(authCtx(t, testAuthToken), &corev1.DeliverRunStartRequest{})
	if status.Code(err) != codes.InvalidArgument {
		t.Fatalf("DeliverRunStart without start_record: err = %v, want InvalidArgument", err)
	}
}

func TestHandleValidation(t *testing.T) {
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	spinup(t, env)
	ctx := authCtx(t, testAuthToken)

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
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	spinup(t, env)
	handle := startRun(t, env)
	ctx := authCtx(t, testAuthToken)

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

	// confirm：释放会话，但绝不关闭 server
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

func TestTwoRunHandleIsolation(t *testing.T) {
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	spinup(t, env)
	handleA, handleB := startRun(t, env), startRun(t, env)
	if handleA == handleB {
		t.Fatal("two runs must receive distinct run_handles")
	}
	ctx := authCtx(t, testAuthToken)

	// 只结束 A：A 报 FINISHED，B 仍是 RUNNING，互不串扰
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

	// confirm A：B 仍可正常走完生命周期
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
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	spinup(t, env)

	const clients = 4
	var wg sync.WaitGroup
	errs := make(chan error, clients)
	for i := 0; i < clients; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			ctx := authCtx(t, testAuthToken)
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
