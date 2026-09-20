package server

import (
	"context"
	"net"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/grpc/metadata"
	"google.golang.org/grpc/status"
	"google.golang.org/grpc/test/bufconn"

	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
)

const (
	bufconnSize      = 1 << 20
	testGrace        = 2 * time.Second
	noShutdownWindow = 200 * time.Millisecond
	callTimeout      = 2 * time.Second

	testOwnerToken = "owner-secret"
	testAuthToken  = "auth-secret"
)

type testEnv struct {
	client corev1.CoreServiceClient
	ctrl   *Controller
}

// newTestEnv 在内存连接上启动完整服务端（含 auth interceptor），返回客户端句柄与关闭控制器。
func newTestEnv(t *testing.T, ownerToken, authToken string) *testEnv {
	t.Helper()
	g := grpc.NewServer(grpc.ChainUnaryInterceptor(UnaryAuthInterceptor(authToken)))
	ctrl := NewController(g, testGrace)
	NewService(ownerToken, authToken, ctrl).Register(g)
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

// authCtx 返回携带 auth token metadata 的带超时 context。
func authCtx(t *testing.T, token string) context.Context {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	t.Cleanup(cancel)
	return metadata.AppendToOutgoingContext(ctx, AuthTokenMetadataKey, token)
}

// plainCtx 返回不带 auth metadata 的带超时 context。
func plainCtx(t *testing.T) context.Context {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	t.Cleanup(cancel)
	return ctx
}

func TestTeardownServiceRejectsWrongToken(t *testing.T) {
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	_, err := env.client.TeardownService(authCtx(t, testAuthToken), &corev1.TeardownServiceRequest{OwnerToken: "wrong-token"})
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
	env := newTestEnv(t, testOwnerToken, testAuthToken)
	if _, err := env.client.TeardownService(authCtx(t, testAuthToken), &corev1.TeardownServiceRequest{OwnerToken: testOwnerToken}); err != nil {
		t.Fatalf("TeardownService: %v", err)
	}
	select {
	case <-env.ctrl.Done():
	case <-time.After(callTimeout):
		t.Fatal("shutdown not completed after valid teardown")
	}
}

func TestTeardownServiceRejectsEmptyConfiguredToken(t *testing.T) {
	env := newTestEnv(t, "", testAuthToken)
	_, err := env.client.TeardownService(authCtx(t, testAuthToken), &corev1.TeardownServiceRequest{OwnerToken: ""})
	if status.Code(err) != codes.PermissionDenied {
		t.Fatalf("TeardownService err = %v, want PermissionDenied", err)
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
	if st := ctrl.Lifecycle().Get(); st != StateClosed {
		t.Fatalf("lifecycle state after shutdown = %v, want StateClosed", st)
	}
}
