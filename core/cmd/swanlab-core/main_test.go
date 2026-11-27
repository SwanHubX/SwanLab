//go:build !windows

package main

import (
	"context"
	"errors"
	"net"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"go.uber.org/goleak"
	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/grpc/status"
	"google.golang.org/grpc/test/bufconn"
	"google.golang.org/protobuf/types/known/wrapperspb"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/portinfo"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	settingsv1 "github.com/swanhubx/swanlab/core/proto/swanlab/settings/core/v1"
)

// staleSocket 留下一个「进程已退出、socket 文件仍在」的 UDS 残留。
func staleSocket(t *testing.T, path string) {
	t.Helper()
	ln, err := net.Listen("unix", path)
	if err != nil {
		t.Fatalf("listen %s: %v", path, err)
	}
	ln.(*net.UnixListener).SetUnlinkOnClose(false)
	if err := ln.Close(); err != nil {
		t.Fatalf("close %s: %v", path, err)
	}
}

// shortDir 返回短路径目录：UDS 路径受 sun_path 限制，t.TempDir 在 macOS 上可能过长。
func shortDir(t *testing.T) string {
	t.Helper()
	dir, err := os.MkdirTemp("/tmp", "slc")
	if err != nil {
		t.Fatalf("mkdtemp: %v", err)
	}
	t.Cleanup(func() { _ = os.RemoveAll(dir) })
	return dir
}

func TestListenAutoUnixReclaimsStaleSocket(t *testing.T) {
	path := filepath.Join(shortDir(t), coreSocketName)
	staleSocket(t, path)

	ln, err := listenAutoUnix(path)
	if err != nil {
		t.Fatalf("listenAutoUnix: %v", err)
	}
	t.Cleanup(func() { _ = ln.Close() })
}

func TestListenAutoUnixRejectsLiveListener(t *testing.T) {
	path := filepath.Join(shortDir(t), coreSocketName)
	live, err := net.Listen("unix", path)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { _ = live.Close() })

	if _, probeErr := listenAutoUnix(path); !errors.Is(probeErr, errEndpointBusy) {
		t.Fatalf("err = %v, want errEndpointBusy", probeErr)
	}
	conn, err := net.Dial("unix", path)
	if err != nil {
		t.Fatalf("live listener broken after probe: %v", err)
	}
	_ = conn.Close()
}

func TestListenAutoUnixKeepsForeignFile(t *testing.T) {
	path := filepath.Join(shortDir(t), coreSocketName)
	content := []byte("caller-owned")
	if err := os.WriteFile(path, content, 0o600); err != nil {
		t.Fatal(err)
	}
	if _, probeErr := listenAutoUnix(path); probeErr == nil || errors.Is(probeErr, errEndpointBusy) {
		t.Fatalf("err = %v, want non-busy error", probeErr)
	}
	got, err := os.ReadFile(path)
	if err != nil || string(got) != string(content) {
		t.Fatalf("foreign file modified: %q err=%v", got, err)
	}
}

func TestOpenEndpointFallsBackOnForeignFile(t *testing.T) {
	dir := shortDir(t)
	if err := os.WriteFile(filepath.Join(dir, coreSocketName), []byte("x"), 0o600); err != nil {
		t.Fatal(err)
	}
	ln, err := openEndpoint("", filepath.Join(dir, "core.port"))
	if err != nil {
		t.Fatalf("openEndpoint: %v", err)
	}
	defer func() { _ = ln.Close() }()
	if _, ok := ln.Addr().(*net.TCPAddr); !ok {
		t.Fatalf("want loopback tcp fallback, got %T", ln.Addr())
	}
}

func TestOpenEndpointBusyWhenUnixLive(t *testing.T) {
	dir := shortDir(t)
	live, err := net.Listen("unix", filepath.Join(dir, coreSocketName))
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { _ = live.Close() })

	if _, probeErr := openEndpoint("", filepath.Join(dir, "core.port")); !errors.Is(probeErr, errEndpointBusy) {
		t.Fatalf("err = %v, want errEndpointBusy", probeErr)
	}
}

func TestOpenEndpointRejectsRelativeUnix(t *testing.T) {
	if _, err := openEndpoint("unix://relative.sock", ""); err == nil {
		t.Fatal("relative unix path must be rejected")
	}
}

func TestOpenEndpointRejectsNonLoopbackTCP(t *testing.T) {
	if _, err := openEndpoint("tcp://0.0.0.0:0", ""); err == nil {
		t.Fatal("non-loopback tcp address must be rejected")
	}
}

func TestPortFileServed(t *testing.T) {
	dir := shortDir(t)
	path := filepath.Join(dir, "core.port")

	if _, alive := portFileServed(path); alive {
		t.Fatal("missing port-file must be treated as free")
	}

	sock := filepath.Join(dir, coreSocketName)
	staleSocket(t, sock)
	info := portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: os.Getpid(), UnixPath: sock}
	if err := portinfo.WriteFile(path, &info); err != nil {
		t.Fatal(err)
	}
	if pid, alive := portFileServed(path); alive || pid != os.Getpid() {
		t.Fatalf("stale endpoint: pid=%d alive=%v, want pid=%d alive=false", pid, alive, os.Getpid())
	}

	// 端点重新可连后，同一份 port-file 应判定为被占用
	if err := os.Remove(sock); err != nil {
		t.Fatal(err)
	}
	live, err := net.Listen("unix", sock)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { _ = live.Close() })
	if _, alive := portFileServed(path); !alive {
		t.Fatal("live endpoint must be reported as served")
	}
}

// TestAssembleServerRegistersRunAndSyncServices 守护宿主接线：Run 与 Sync 两个 gRPC service
// 都必须在 assembleServer 里注册，删掉任意一行都会让 CLI 侧 RPC 直接 NOT_FOUND。
func TestAssembleServerRegistersRunAndSyncServices(t *testing.T) {
	grpcServer, _, mgr, syncMgr, _ := assembleServer(0)
	t.Cleanup(func() {
		grpcServer.Stop()
		_ = mgr.Close()
		_ = syncMgr.Close()
	})

	info := grpcServer.GetServiceInfo()
	for _, name := range []string{
		corev1.CoreService_ServiceDesc.ServiceName,
		corev1.CoreSyncService_ServiceDesc.ServiceName,
	} {
		if _, ok := info[name]; !ok {
			t.Fatalf("gRPC service %q is not registered; got %v", name, info)
		}
	}
}

type panicTestService interface {
	Call(context.Context, *wrapperspb.StringValue) (*wrapperspb.StringValue, error)
}

type panicTestHandler struct{}

func (panicTestHandler) Call(_ context.Context, req *wrapperspb.StringValue) (*wrapperspb.StringValue, error) {
	if req.GetValue() == "panic" {
		panic("boom")
	}
	return req, nil
}

var panicTestServiceDesc = grpc.ServiceDesc{
	ServiceName: "test.PanicService",
	HandlerType: (*panicTestService)(nil),
	Methods: []grpc.MethodDesc{{
		MethodName: "Call",
		Handler: func(srv any, ctx context.Context, dec func(any) error, interceptor grpc.UnaryServerInterceptor) (any, error) {
			in := new(wrapperspb.StringValue)
			if err := dec(in); err != nil {
				return nil, err
			}
			if interceptor == nil {
				return srv.(panicTestService).Call(ctx, in)
			}
			info := &grpc.UnaryServerInfo{Server: srv, FullMethod: "/test.PanicService/Call"}
			return interceptor(ctx, in, info, func(ctx context.Context, req any) (any, error) {
				return srv.(panicTestService).Call(ctx, req.(*wrapperspb.StringValue))
			})
		},
	}},
}

func TestAssembleServerRecoversHandlerPanic(t *testing.T) {
	grpcServer, _, mgr, syncMgr, _ := assembleServer(0)
	grpcServer.RegisterService(&panicTestServiceDesc, panicTestHandler{})
	lis := bufconn.Listen(1 << 20)
	t.Cleanup(func() {
		grpcServer.Stop()
		_ = lis.Close()
		_ = mgr.Close()
		_ = syncMgr.Close()
	})
	go func() { _ = grpcServer.Serve(lis) }()

	conn, err := grpc.NewClient(
		"passthrough:///bufnet",
		grpc.WithContextDialer(func(ctx context.Context, _ string) (net.Conn, error) {
			return lis.DialContext(ctx)
		}),
		grpc.WithTransportCredentials(insecure.NewCredentials()),
	)
	if err != nil {
		t.Fatalf("new client: %v", err)
	}
	t.Cleanup(func() { _ = conn.Close() })

	ctx, cancel := context.WithTimeout(context.Background(), 2*time.Second)
	defer cancel()
	out := new(wrapperspb.StringValue)
	err = conn.Invoke(ctx, "/test.PanicService/Call", wrapperspb.String("panic"), out)
	if status.Code(err) != codes.Internal {
		t.Fatalf("panic err=%v, want codes.Internal", err)
	}
	if !strings.Contains(err.Error(), "internal panic: boom") {
		t.Fatalf("panic err=%v, want internal panic message", err)
	}

	err = conn.Invoke(ctx, "/test.PanicService/Call", wrapperspb.String("ok"), out)
	if err != nil || out.GetValue() != "ok" {
		t.Fatalf("normal call resp=%v err=%v, want ok/nil", out, err)
	}
}

type sessionPanicTestService interface {
	Call(context.Context, *corev1.GetOperationStatsRequest) (*wrapperspb.StringValue, error)
}

type sessionPanicTestHandler struct{}

func (sessionPanicTestHandler) Call(context.Context, *corev1.GetOperationStatsRequest) (*wrapperspb.StringValue, error) {
	panic("session boom")
}

var sessionPanicTestServiceDesc = grpc.ServiceDesc{
	ServiceName: "test.SessionPanicService",
	HandlerType: (*sessionPanicTestService)(nil),
	Methods: []grpc.MethodDesc{{
		MethodName: "Call",
		Handler: func(srv any, ctx context.Context, dec func(any) error, interceptor grpc.UnaryServerInterceptor) (any, error) {
			in := new(corev1.GetOperationStatsRequest)
			if err := dec(in); err != nil {
				return nil, err
			}
			if interceptor == nil {
				return srv.(sessionPanicTestService).Call(ctx, in)
			}
			info := &grpc.UnaryServerInfo{Server: srv, FullMethod: "/test.SessionPanicService/Call"}
			return interceptor(ctx, in, info, func(ctx context.Context, req any) (any, error) {
				return srv.(sessionPanicTestService).Call(ctx, req.(*corev1.GetOperationStatsRequest))
			})
		},
	}},
}

func TestAssembleServerPanicInvalidatesActiveSession(t *testing.T) {
	baseline := goleak.IgnoreCurrent()
	grpcServer, ctrl, mgr, syncMgr, _ := assembleServer(0)
	grpcServer.RegisterService(&sessionPanicTestServiceDesc, sessionPanicTestHandler{})
	lis := bufconn.Listen(1 << 20)
	serveDone := make(chan struct{})
	go func() {
		defer close(serveDone)
		_ = grpcServer.Serve(lis)
	}()

	conn, err := grpc.NewClient(
		"passthrough:///bufnet",
		grpc.WithContextDialer(func(ctx context.Context, _ string) (net.Conn, error) {
			return lis.DialContext(ctx)
		}),
		grpc.WithTransportCredentials(insecure.NewCredentials()),
	)
	if err != nil {
		t.Fatalf("new client: %v", err)
	}

	ctx, cancel := context.WithTimeout(context.Background(), 5*time.Second)
	defer cancel()
	coreClient := corev1.NewCoreServiceClient(conn)
	if _, err = coreClient.SpinupService(ctx, &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
	startResp, err := coreClient.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: &settingsv1.CoreSettings{
			RunId:          "panic-run",
			RunDir:         t.TempDir(),
			Mode:           settingsv1.CoreMode_CORE_MODE_LOCAL,
			RecordBatch:    100,
			RecordInterval: 1.5,
			SaveSize:       50 << 30,
			SaveSplit:      100 << 20,
			SavePart:       32 << 20,
			SaveBatch:      100,
		},
		StartRecord: &runv1.StartRecord{Project: "demo", Id: "panic-run"},
	})
	if err != nil || !startResp.GetSuccess() || startResp.GetSessionId() == "" {
		t.Fatalf("DeliverRunStart resp=%v err=%v", startResp, err)
	}
	sessionID := startResp.GetSessionId()

	out := new(wrapperspb.StringValue)
	err = conn.Invoke(ctx, "/test.SessionPanicService/Call", &corev1.GetOperationStatsRequest{SessionId: sessionID}, out)
	if status.Code(err) != codes.Internal {
		t.Fatalf("panic err=%v, want codes.Internal", err)
	}

	if _, err = coreClient.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: sessionID}); status.Code(err) != codes.NotFound {
		t.Fatalf("GetOperationStats after panic err=%v, want codes.NotFound", err)
	}

	_ = conn.Close()
	ctrl.Shutdown("test")
	<-ctrl.Done()
	<-serveDone
	_ = lis.Close()
	_ = mgr.Close()
	_ = syncMgr.Close()
	goleak.VerifyNone(t, baseline)
}

func TestPanicRecoveryInterceptorDeduplicatesStackByMethod(t *testing.T) {
	logDir := t.TempDir()
	console.Reset()
	if err := console.Init(logDir); err != nil {
		t.Fatalf("console.Init: %v", err)
	}
	t.Cleanup(console.Reset)

	var invalidated []string
	interceptor := newPanicRecoveryInterceptor(func(id string) {
		invalidated = append(invalidated, id)
	})

	call := func(method string, req any) error {
		_, err := interceptor(context.Background(), req, &grpc.UnaryServerInfo{FullMethod: method}, func(context.Context, any) (any, error) {
			panic("boom")
		})
		return err
	}

	req := &corev1.GetOperationStatsRequest{SessionId: "sess-1"}
	if err := call("/svc/MethodA", req); status.Code(err) != codes.Internal {
		t.Fatalf("MethodA #1 err=%v, want Internal", err)
	}
	if err := call("/svc/MethodA", req); status.Code(err) != codes.Internal {
		t.Fatalf("MethodA #2 err=%v, want Internal", err)
	}
	if err := call("/svc/MethodB", wrapperspb.String("no-session")); status.Code(err) != codes.Internal {
		t.Fatalf("MethodB #1 err=%v, want Internal", err)
	}

	if len(invalidated) != 2 || invalidated[0] != "sess-1" || invalidated[1] != "sess-1" {
		t.Fatalf("invalidated = %v, want [sess-1 sess-1]", invalidated)
	}

	data, err := os.ReadFile(filepath.Join(logDir, "debug-core.log"))
	if err != nil {
		t.Fatalf("read debug-core.log: %v", err)
	}
	text := string(data)
	if got := strings.Count(text, "stack suppressed"); got != 1 {
		t.Fatalf("stack suppressed count = %d, want 1\nlog:\n%s", got, text)
	}
	if got := strings.Count(text, "runtime/debug.Stack()"); got != 2 {
		t.Fatalf("full stack trace count = %d, want 2 (once per method)\nlog:\n%s", got, text)
	}
}
