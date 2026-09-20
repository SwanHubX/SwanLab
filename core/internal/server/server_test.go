package server

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

	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
)

const (
	bufconnSize      = 1 << 20
	testGrace        = 2 * time.Second
	noShutdownWindow = 200 * time.Millisecond
	callTimeout      = 2 * time.Second
)

type testEnv struct {
	client corev1.CoreServiceClient
	ctrl   *Controller
}

// newTestEnv 在内存连接上启动完整服务端，返回客户端句柄与关闭控制器。
func newTestEnv(t *testing.T, ownerToken string) *testEnv {
	t.Helper()
	g := grpc.NewServer()
	ctrl := NewController(g, testGrace)
	NewService(ownerToken, ctrl).Register(g)
	lis := bufconn.Listen(bufconnSize)
	go func() { _ = g.Serve(lis) }()
	conn, err := grpc.NewClient("passthrough:///bufnet",
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

func TestTeardownServiceRejectsWrongToken(t *testing.T) {
	env := newTestEnv(t, "owner-secret")
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	defer cancel()
	_, err := env.client.TeardownService(ctx, &corev1.TeardownServiceRequest{OwnerToken: "wrong-token"})
	if status.Code(err) != codes.PermissionDenied {
		t.Fatalf("TeardownService err = %v, want PermissionDenied", err)
	}
	select {
	case <-env.ctrl.Done():
		t.Fatal("wrong owner token must not trigger shutdown")
	case <-time.After(noShutdownWindow):
	}
}

func TestTeardownServiceShutsDownServer(t *testing.T) {
	env := newTestEnv(t, "owner-secret")
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	defer cancel()
	if _, err := env.client.TeardownService(ctx, &corev1.TeardownServiceRequest{OwnerToken: "owner-secret"}); err != nil {
		t.Fatalf("TeardownService: %v", err)
	}
	select {
	case <-env.ctrl.Done():
	case <-time.After(callTimeout):
		t.Fatal("shutdown not completed after valid teardown")
	}
}

func TestTeardownServiceRejectsEmptyConfiguredToken(t *testing.T) {
	env := newTestEnv(t, "")
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	defer cancel()
	_, err := env.client.TeardownService(ctx, &corev1.TeardownServiceRequest{OwnerToken: ""})
	if status.Code(err) != codes.PermissionDenied {
		t.Fatalf("TeardownService err = %v, want PermissionDenied", err)
	}
}

func TestRunLevelRPCsUnimplemented(t *testing.T) {
	env := newTestEnv(t, "owner-secret")
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	defer cancel()
	if _, err := env.client.GetCapabilities(ctx, &corev1.GetCapabilitiesRequest{}); status.Code(err) != codes.Unimplemented {
		t.Fatalf("GetCapabilities err = %v, want Unimplemented", err)
	}
	if _, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{}); status.Code(err) != codes.Unimplemented {
		t.Fatalf("UpsertScalars err = %v, want Unimplemented", err)
	}
	if _, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{}); status.Code(err) != codes.Unimplemented {
		t.Fatalf("DeliverRunStart err = %v, want Unimplemented", err)
	}
}

func TestControllerShutdownIdempotent(t *testing.T) {
	g := grpc.NewServer()
	ctrl := NewController(g, testGrace)
	ctrl.Shutdown("first")
	ctrl.Shutdown("second")
	select {
	case <-ctrl.Done():
	case <-time.After(callTimeout):
		t.Fatal("controller Done not closed after Shutdown")
	}
	// 重复读取已关闭的 Done 不应阻塞或 panic
	<-ctrl.Done()
}
