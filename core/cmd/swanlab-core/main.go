// Command swanlab-core 是 SwanLab Go core 的进程入口。
//
// 端点约定：
//
//	--listen unix:///path/to/uds   Linux/macOS 进程内通信（手动调试入口）
//	--listen tcp://127.0.0.1:port  Windows 回环地址
//	--port-filename <路径>         listen 成功后原子写入端点回报文件（SDK 启动约定）
//
// 未传 --listen 时按平台自选端点：POSIX 使用 port-filename 同目录下的
// core.sock（UDS，目录需已存在），Windows 使用 127.0.0.1 随机回环端口。
// --port-filename 与 --owner-token-file 成对出现，owner token 是唯一允许
// 触发服务级关闭的凭证，通过私有文件传入，不得出现在命令行或日志中。
//
// 生命周期：Teardown RPC、SIGINT/SIGTERM、父进程退出（process 包监控）或
// Serve 异常统一汇入 service controller 的关闭路径（GracefulStop → 超时
// 强制 Stop）；退出时只清理自己创建的 socket 文件与 port-file。
package main

import (
	"context"
	"errors"
	"flag"
	"fmt"
	"net"
	"os"
	"os/signal"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"syscall"
	"time"

	"google.golang.org/grpc"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/portinfo"
	"github.com/swanhubx/swanlab/core/internal/pkg/process"
	"github.com/swanhubx/swanlab/core/internal/server"
)

// version 与 commit 由构建管线通过 -ldflags -X 注入（见 core/hatch.py），
// 缺省值仅供本地 go run / go build 使用。
var (
	version = "dev"
	commit  = "unknown"
)

// 与命令行参数等价的环境变量，供 Python SDK spawn 时注入。
const (
	envListenAddr = "SWANLAB_CORE_LISTEN"
	envParentPID  = "SWANLAB_CORE_PARENT_PID"
)

// 进程退出码约定：0 正常退出（含信号触发的优雅关闭）；2 用法错误；1 运行错误。
const (
	exitUsageError = 2
	exitRunError   = 1
)

// 自选端点与收尾参数。
const (
	coreSocketName    = "core.sock"
	loopbackAddr      = "127.0.0.1:0"
	shutdownGrace     = 10 * time.Second
	secretFileMaxSize = 4096
)

func main() {
	os.Exit(run(os.Args[1:]))
}

func run(args []string) int {
	fs := flag.NewFlagSet("swanlab-core", flag.ContinueOnError)
	printVersion := fs.Bool("version", false, "print version and exit")
	listenAddr := fs.String("listen", os.Getenv(envListenAddr),
		"listen endpoint, unix://<uds path> or tcp://<addr:port>; auto-selected per platform when unset")
	portFilename := fs.String("port-filename", "",
		"endpoint report file; atomically written once listen succeeds, for callers to poll")
	ownerTokenFile := fs.String("owner-token-file", "",
		"owner token file, held only by the service owner; sole credential for service-level teardown")
	parentPID := fs.Int("parent-pid", envInt(envParentPID),
		"expected parent PID; core exits when the parent exits, defaults to the actual parent at startup")
	if err := fs.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			return 0
		}
		return exitUsageError
	}

	if *printVersion {
		fmt.Printf("swanlab-core %s (commit %s)\n", version, commit)
		return 0
	}
	if *listenAddr == "" && *portFilename == "" {
		console.Error("no listen endpoint: pass --listen (manual debug) or --port-filename (SDK startup convention)")
		return exitUsageError
	}
	if (*portFilename != "") != (*ownerTokenFile != "") {
		console.Error("--port-filename and --owner-token-file must be provided together")
		return exitUsageError
	}

	// 自建资源记录，退出时只清理自己创建的部分。
	var socketPath string
	wrotePortFile := false
	defer func() {
		cleanupSocket(socketPath)
		if wrotePortFile {
			_ = os.Remove(*portFilename)
		}
	}()

	// owner token 先于任何资源创建读取，尽早失败。
	ownerToken, err := readOwnerToken(*ownerTokenFile)
	if err != nil {
		console.Error("failed to read owner token:", err)
		return exitRunError
	}

	ln, err := openEndpoint(*listenAddr, *portFilename)
	if err != nil {
		console.Error("listen failed:", err)
		return exitRunError
	}
	defer func() { _ = ln.Close() }()
	if addr, ok := ln.Addr().(*net.UnixAddr); ok && !strings.HasPrefix(addr.Name, "@") {
		socketPath = addr.Name
	}

	// 父进程监控：显式传入的 PID 优先（启动约定），未传时回退为监控启动
	// 瞬间的实际父进程（本地终端运行场景）。监控建立失败按约定终止启动。
	pid := *parentPID
	if pid <= 0 {
		pid = os.Getppid()
	}
	parentExited, err := process.NotifyOnParentExit(pid)
	if err != nil {
		console.Error("failed to watch parent process, aborting startup:", err)
		return exitRunError
	}

	grpcServer := grpc.NewServer()
	ctrl := server.NewController(grpcServer, shutdownGrace)
	server.NewService(ownerToken, ctrl).Register(grpcServer)

	// listen 与 server 初始化均成功后才写 port-file。
	if *portFilename != "" {
		authToken, err2 := portinfo.NewAuthToken()
		if err2 != nil {
			console.Error("failed to generate auth token:", err2)
			return exitRunError
		}
		info := portinfo.Info{Protocol: portinfo.ProtocolVersion, AuthToken: authToken}
		switch addr := ln.Addr().(type) {
		case *net.UnixAddr:
			info.UnixPath = addr.Name
		case *net.TCPAddr:
			info.SockPort = addr.Port
		default:
			console.Error("unrecognized listener address type:", ln.Addr())
			return exitRunError
		}
		if err2 = portinfo.WriteFile(*portFilename, &info); err2 != nil {
			console.Error("failed to write port-file:", err2)
			return exitRunError
		}
		wrotePortFile = true
	}

	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer stop()

	console.Infof("swanlab-core %s listening on %s (parent pid %d)", version, ln.Addr(), pid)

	serveErr := make(chan error, 1)
	go func() {
		serveErr <- grpcServer.Serve(ln)
	}()

	exitCode := 0
	var cause string
	var serveFailure error
	select {
	case <-ctx.Done():
		console.Info("shutdown signal received, stopping")
		cause = "signal"
	case <-parentExited:
		console.Warning("parent process exited, stopping core")
		cause = "parent-exit"
	case err := <-serveErr:
		serveFailure = err
		cause = "serve-error"
	}
	ctrl.Shutdown(cause)
	// 等待 Serve 返回与关闭序列完成（两者任一先行均可）。
	select {
	case serveFailure = <-serveErr:
	case <-ctrl.Done():
	}
	<-ctrl.Done()
	if serveFailure != nil {
		console.Error("gRPC Serve exited with error:", serveFailure)
		exitCode = exitRunError
	}
	return exitCode
}

// openEndpoint 创建监听器。显式 --listen 优先（手动调试）；否则按平台自选：
// POSIX 使用 port-filename 同目录下的 UDS，Windows 使用随机回环端口。
func openEndpoint(listenAddr, portFilename string) (net.Listener, error) {
	if listenAddr == "" {
		if runtime.GOOS == "windows" {
			return net.Listen("tcp", loopbackAddr)
		}
		sockPath := filepath.Join(filepath.Dir(portFilename), coreSocketName)
		return net.Listen("unix", sockPath)
	}
	scheme, rest, ok := strings.Cut(listenAddr, "://")
	if !ok {
		return nil, fmt.Errorf("listen endpoint missing scheme prefix (unix:// or tcp://): %s", listenAddr)
	}
	switch scheme {
	case "unix":
		if runtime.GOOS == "windows" {
			return nil, errors.New("unix:// endpoints are not supported on Windows; use tcp://127.0.0.1:<port>")
		}
		return net.Listen("unix", rest)
	case "tcp":
		return net.Listen("tcp", rest)
	default:
		return nil, fmt.Errorf("unsupported listen scheme %q (only unix:// or tcp://)", scheme)
	}
}

// readOwnerToken 读取 owner token 文件；路径为空返回空串。内容去除首尾空白，
// 拒绝空文件与超长文件，token 值不进入日志。
func readOwnerToken(path string) (string, error) {
	if path == "" {
		return "", nil
	}
	data, err := os.ReadFile(path)
	if err != nil {
		return "", fmt.Errorf("read %s: %w", path, err)
	}
	if len(data) > secretFileMaxSize {
		return "", fmt.Errorf("owner token file exceeds %d bytes: %s", secretFileMaxSize, path)
	}
	token := strings.TrimSpace(string(data))
	if token == "" {
		return "", fmt.Errorf("owner token file is empty: %s", path)
	}
	return token, nil
}

// cleanupSocket 删除自己创建的 UDS socket 文件；Linux 抽象 socket（@ 前缀）
// 不占文件系统，无需清理。
func cleanupSocket(path string) {
	if path == "" || strings.HasPrefix(path, "@") {
		return
	}
	_ = os.Remove(path)
}

// envInt 解析整型环境变量，缺失或非法时返回 0。
func envInt(name string) int {
	v, err := strconv.Atoi(strings.TrimSpace(os.Getenv(name)))
	if err != nil {
		return 0
	}
	return v
}
