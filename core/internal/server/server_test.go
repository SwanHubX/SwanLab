package server

import (
	"context"
	"net"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/credentials/insecure"
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
func newTestEnv(t *testing.T) *testEnv {
	t.Helper()
	g := grpc.NewServer()
	ctrl := NewController(g, testGrace)
	NewService(ctrl).Register(g)
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

// callCtx 返回带超时的 context。
func callCtx(t *testing.T) context.Context {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	t.Cleanup(cancel)
	return ctx
}

func TestTeardownServiceShutsDownServer(t *testing.T) {
	env := newTestEnv(t)
	if _, err := env.client.TeardownService(callCtx(t), &corev1.TeardownServiceRequest{}); err != nil {
		t.Fatalf("TeardownService: %v", err)
	}
	select {
	case <-env.ctrl.Done():
	case <-time.After(callTimeout):
		t.Fatal("shutdown not completed after valid teardown")
	}
}

// TestTeardownServiceIgnoresOwnerToken 不校验 owner_token：proto 字段保留但忽略。
func TestTeardownServiceIgnoresOwnerToken(t *testing.T) {
	env := newTestEnv(t)
	if _, err := env.client.TeardownService(callCtx(t), &corev1.TeardownServiceRequest{OwnerToken: "ignored"}); err != nil {
		t.Fatalf("TeardownService: %v", err)
	}
	select {
	case <-env.ctrl.Done():
	case <-time.After(callTimeout):
		t.Fatal("shutdown not completed after teardown")
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
