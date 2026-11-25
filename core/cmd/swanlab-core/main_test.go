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

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/grpc/status"
	"google.golang.org/grpc/test/bufconn"
	"google.golang.org/protobuf/types/known/wrapperspb"

	"github.com/swanhubx/swanlab/core/internal/pkg/portinfo"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
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
