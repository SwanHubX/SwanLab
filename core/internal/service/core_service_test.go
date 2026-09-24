package service

import (
	"context"
	"net"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/grpc/status"
	"google.golang.org/grpc/test/bufconn"

	"github.com/swanhubx/swanlab/core/internal/manager"
	"github.com/swanhubx/swanlab/core/internal/server"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

const (
	bufconnSize      = 1 << 20
	callTimeout      = 2 * time.Second
	noShutdownWindow = 200 * time.Millisecond
)

type testEnv struct {
	client corev1.CoreServiceClient
	ctrl   *server.Controller
}

func newTestEnv(t *testing.T) *testEnv {
	t.Helper()
	g := grpc.NewServer()
	ctrl := server.NewController(g, 2*time.Second)
	NewCoreService(ctrl, manager.New()).Register(g)
	lis := bufconn.Listen(bufconnSize)
	go func() { _ = g.Serve(lis) }()
	conn, err := grpc.NewClient(
		"passthrough:///bufnet",
		grpc.WithTransportCredentials(insecure.NewCredentials()),
		grpc.WithContextDialer(func(context.Context, string) (net.Conn, error) {
			return lis.DialContext(context.Background())
		}),
	)
	if err != nil {
		t.Fatalf("grpc.NewClient: %v", err)
	}
	t.Cleanup(func() {
		_ = conn.Close()
		ctrl.Shutdown("test-cleanup")
		<-ctrl.Done()
	})
	return &testEnv{client: corev1.NewCoreServiceClient(conn), ctrl: ctrl}
}

func callCtx(t *testing.T) context.Context {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	t.Cleanup(cancel)
	return ctx
}

func spinup(t *testing.T, env *testEnv) {
	t.Helper()
	if _, err := env.client.SpinupService(callCtx(t), &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
}

func startRun(t *testing.T, env *testEnv) string {
	t.Helper()
	resp, err := env.client.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{
		StartRecord: &runv1.StartRecord{},
	})
	if err != nil || !resp.GetSuccess() || resp.GetRunHandle() == "" {
		t.Fatalf("DeliverRunStart: resp=%v err=%v", resp, err)
	}
	return resp.GetRunHandle()
}

func TestSpinupAndTeardown(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	spinup(t, env)
	if _, err := env.client.TeardownService(callCtx(t), &corev1.TeardownServiceRequest{}); err != nil {
		t.Fatalf("TeardownService: %v", err)
	}
	select {
	case <-env.ctrl.Done():
	case <-time.After(callTimeout):
		t.Fatal("shutdown not completed after teardown")
	}
}

func TestRunRPCsRejectedBeforeReady(t *testing.T) {
	env := newTestEnv(t)
	ctx := callCtx(t)
	for name, call := range map[string]func() error{
		"start": func() error {
			_, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{StartRecord: &runv1.StartRecord{}})
			return err
		},
		"upsert": func() error {
			_, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{RunHandle: "any"})
			return err
		},
		"stats": func() error {
			_, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: "any"})
			return err
		},
		"finish": func() error {
			_, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: "any"})
			return err
		},
		"confirm": func() error {
			_, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: "any"})
			return err
		},
	} {
		if err := call(); status.Code(err) != codes.FailedPrecondition {
			t.Fatalf("%s before ready: err=%v, want FailedPrecondition", name, err)
		}
	}
}

func TestUpsertsRemainUnimplemented(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env)
	ctx := callCtx(t)
	calls := []func() error{
		func() error {
			_, err := env.client.UpsertColumns(ctx, &corev1.UpsertColumnsRequest{RunHandle: handle})
			return err
		},
		func() error {
			_, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{RunHandle: handle})
			return err
		},
		func() error {
			_, err := env.client.UpsertMedia(ctx, &corev1.UpsertMediaRequest{RunHandle: handle})
			return err
		},
		func() error {
			_, err := env.client.UpsertLogs(ctx, &corev1.UpsertLogsRequest{RunHandle: handle})
			return err
		},
		func() error {
			_, err := env.client.UpsertSaves(ctx, &corev1.UpsertSavesRequest{RunHandle: handle})
			return err
		},
	}
	for i, call := range calls {
		if err := call(); status.Code(err) != codes.Unimplemented {
			t.Fatalf("upsert %d err=%v, want Unimplemented", i, err)
		}
	}
}

func TestRunLifecycleAndErrorMapping(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)

	if _, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("missing start_record err=%v, want InvalidArgument", err)
	}
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("empty handle err=%v, want InvalidArgument", err)
	}
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: "missing"}); status.Code(err) != codes.NotFound {
		t.Fatalf("unknown handle err=%v, want NotFound", err)
	}

	handle := startRun(t, env)
	resp, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle})
	if err != nil || resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("running stats: resp=%v err=%v", resp, err)
	}
	if _, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); status.Code(err) != codes.FailedPrecondition {
		t.Fatalf("early confirm err=%v, want FailedPrecondition", err)
	}
	if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{RunHandle: handle}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	resp, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle})
	if err != nil || resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_FINISHED {
		t.Fatalf("finished stats: resp=%v err=%v", resp, err)
	}
	if _, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{RunHandle: handle}); err != nil {
		t.Fatalf("ConfirmRunFinish: %v", err)
	}
	select {
	case <-env.ctrl.Done():
		t.Fatal("ConfirmRunFinish must not stop the server")
	case <-time.After(noShutdownWindow):
	}
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{RunHandle: handle}); status.Code(err) != codes.NotFound {
		t.Fatalf("released handle err=%v, want NotFound", err)
	}
}
