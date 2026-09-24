// Command swanlab-core 是 SwanLab Go core 的进程入口。
//
// 启动约定（Owner Mode）：
//
//	swanlab-core --port-filename <runtime>/core.port --parent-pid <pid>
//
// 端点约定：
//
//	--listen unix:///path/to/uds   Linux/macOS 进程内通信（手动调试入口）
//	--listen tcp://127.0.0.1:port  Windows 回环地址
//	--port-filename <路径>         listen 成功后原子写入端点回报文件（SDK 启动约定）
//
// 未传 --listen 时按平台自选端点：POSIX 使用 port-filename 同目录下的
// core.sock（UDS，目录需已存在），listen 失败记录 warning 后回退
// 127.0.0.1 随机回环端口；Windows 使用随机回环端口。
//
// 信任模型：不使用应用层 token。POSIX 使用 UDS（socket 位于 owner-only
// 私有 runtime 目录），port-file 以 0600 原子发布；loopback TCP 回退
// 不隔离本机用户，能连接端口的本机进程可调用全部 RPC（含 Teardown）。
//
// --detach 与 --idle-timeout 为 detached 模式预留：detached 未实现，
// 传入时忽略（no-op）并记录 warning，core 仍以 owner 模式运行。
//
// 生命周期：Teardown RPC、SIGINT/SIGTERM、父进程退出（process 包监控）或
// Serve 异常汇入 service controller 的关闭路径（GracefulStop → 超时强制
// Stop）；退出时清理自己创建的 socket 文件与 port-file。
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
// 缺省值供本地 go run / go build 使用。
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
	coreSocketName = "core.sock"
	loopbackAddr   = "127.0.0.1:0"
	shutdownGrace  = 10 * time.Second
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
	parentPID := fs.Int("parent-pid", envInt(envParentPID),
		"expected parent PID; core exits when the parent exits, defaults to the actual parent at startup")
	detach := fs.Bool("detach", false,
		"detached mode (reserved; accepted but ignored in this build)")
	idleTimeout := fs.Duration("idle-timeout", 0,
		"detached idle timeout (reserved; accepted but ignored in this build)")
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
	if *detach || *idleTimeout != 0 {
		console.Warning("--detach/--idle-timeout accepted but ignored: detached mode is not implemented; running in owner mode")
	}
	if *listenAddr == "" && *portFilename == "" {
		console.Error("no listen endpoint: pass --listen (manual debug) or --port-filename (SDK startup convention)")
		return exitUsageError
	}

	// 自建资源记录，退出时清理自己创建的部分。
	var socketPath string
	selfPID := os.Getpid() // 退出时凭它确认 port-file 属本实例
	wrotePortFile := false
	defer func() {
		cleanupSocket(socketPath)
		if wrotePortFile {
			cleanupPortFile(*portFilename, selfPID)
		}
	}()

	ln, err := openEndpoint(*listenAddr, *portFilename)
	if err != nil {
		console.Error("listen failed:", err)
		return exitRunError
	}
	defer func() { _ = ln.Close() }()
	if addr, ok := ln.Addr().(*net.UnixAddr); ok && !strings.HasPrefix(addr.Name, "@") {
		socketPath = addr.Name
	}

	// 父进程监控：显式传入的 PID 生效（启动约定）；未传时监控启动瞬间的
	// 实际父进程（本地终端运行场景）。监控建立失败终止启动。
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
	server.NewService(ctrl).Register(grpcServer)

	// listen 与 server 初始化成功后写 port-file。
	if *portFilename != "" {
		info := portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: selfPID}
		switch addr := ln.Addr().(type) {
		case *net.UnixAddr:
			info.UnixPath = addr.Name
		case *net.TCPAddr:
			info.SockPort = addr.Port
		default:
			console.Error("unrecognized listener address type:", ln.Addr())
			return exitRunError
		}
		if err2 := portinfo.WriteFile(*portFilename, &info); err2 != nil {
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
	served := false
	select {
	case <-ctx.Done():
		console.Info("shutdown signal received, stopping")
		cause = "signal"
	case <-parentExited:
		console.Warning("parent process exited, stopping core")
		cause = "parent-exit"
	case serveFailure = <-serveErr:
		served = true
		cause = "serve-error"
	}
	ctrl.Shutdown(cause)
	// Shutdown 的关闭序列会让 Serve 返回（GracefulStop 超时转 Stop）；
	// 上面的 select 未消费 serveErr 时，在此收取其结果。
	if !served {
		serveFailure = <-serveErr
	}
	<-ctrl.Done()
	if serveFailure != nil {
		console.Error("gRPC Serve exited with error:", serveFailure)
		exitCode = exitRunError
	}
	return exitCode
}

// openEndpoint 创建监听器。显式 --listen 生效（手动调试，失败不回退，
// tcp:// 限定回环）；未传时按平台自选：POSIX 尝试 port-filename 同目录
// 下的 UDS，失败记录 warning 后回退随机回环 TCP；Windows 使用随机回环端口。
func openEndpoint(listenAddr, portFilename string) (net.Listener, error) {
	if listenAddr == "" {
		if runtime.GOOS == "windows" {
			return net.Listen("tcp", loopbackAddr)
		}
		sockPath := filepath.Join(filepath.Dir(portFilename), coreSocketName)
		ln, err := net.Listen("unix", sockPath)
		if err != nil {
			console.Warningf("unix listen on %s failed: %v; falling back to loopback tcp %s", sockPath, err, loopbackAddr)
			return net.Listen("tcp", loopbackAddr)
		}
		return ln, nil
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
		if err := requireLoopbackHost(rest); err != nil {
			return nil, err
		}
		return net.Listen("tcp", rest)
	default:
		return nil, fmt.Errorf("unsupported listen scheme %q (only unix:// or tcp://)", scheme)
	}
}

// cleanupSocket 删除自己创建的 UDS socket 文件；Linux 抽象 socket（@ 前缀）
// 不占文件系统，无需清理。
func cleanupSocket(path string) {
	if path == "" || strings.HasPrefix(path, "@") {
		return
	}
	_ = os.Remove(path)
}

// cleanupPortFile 在 pid 匹配时删除 port-file。同一路径被后启动实例
// 覆盖后，本实例退出不得误删他者文件；文件缺失、损坏或 pid 不匹配时不删。
func cleanupPortFile(path string, pid int) {
	if path == "" || pid <= 0 {
		return
	}
	info, err := portinfo.ParseFile(path)
	if err != nil {
		return
	}
	if info.PID != pid {
		return
	}
	_ = os.Remove(path)
}

// requireLoopbackHost 限制显式 tcp:// 监听地址为回环，避免服务暴露到
// 网络；放行 127.0.0.1/::1/localhost。
func requireLoopbackHost(addr string) error {
	host, _, err := net.SplitHostPort(addr)
	if err != nil {
		return fmt.Errorf("parse tcp listen address %q: %w", addr, err)
	}
	if host == "localhost" {
		return nil
	}
	if ip := net.ParseIP(host); ip != nil && ip.IsLoopback() {
		return nil
	}
	return fmt.Errorf("tcp listen address %q is not loopback; use 127.0.0.1 or ::1", addr)
}

// envInt 解析整型环境变量，缺失或非法时返回 0。
func envInt(name string) int {
	v, err := strconv.Atoi(strings.TrimSpace(os.Getenv(name)))
	if err != nil {
		return 0
	}
	return v
}
