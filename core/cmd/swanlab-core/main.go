// Command swanlab-core 是 SwanLab Go core 的进程入口。
//
// 启动与端点约定：
//
//	--port-filename <路径>         listen 成功后原子写入端点回报文件（SDK 发现端点依据）
//	--parent-pid <pid>             监控指定的父进程 PID，父进程退出则 core 联动退出
//	--detach                       脱离父进程运行，直到收到信号、Teardown 或空闲超时
//	--idle-timeout <时长>          活跃连接归零后的空闲退出等待时间（0 表示不启用）
//	--listen unix:///path/to/uds   Linux/macOS 本地通信（手动调试入口）
//	--listen tcp://127.0.0.1:port  Windows 回环地址
//	--log-dir <目录>               诊断日志目录（debug-core.log），为空时只写 stderr
//
// 未显式指定 --listen 时按平台自选端点：
// POSIX 优先使用 port-filename 同目录下的 core.sock（UDS），listen 失败时记录 warning 并回退至 127.0.0.1 随机回环端口；
// Windows 默认直接使用随机回环端口。
//
// 端点占用处理：
// 已有实例在服务时直接失败（退出码 3），不切换传输、不覆盖 port-file；
// 崩溃残留的 socket 回收后重建，其它文件不删除、回退回环 TCP。
//
// 生命周期管理：
// 汇集 Teardown RPC、系统信号、父进程退出监控或 Serve 异常，统一触发 Controller 的
// 关闭序列（优先 GracefulStop，超时强制 Stop）；退出时清理自己创建的 socket 与 port-file 文件。
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
	"runtime/debug"
	"strconv"
	"strings"
	"sync"
	"syscall"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/manager"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/portinfo"
	"github.com/swanhubx/swanlab/core/internal/pkg/process"
	"github.com/swanhubx/swanlab/core/internal/server"
	"github.com/swanhubx/swanlab/core/internal/service"
	"github.com/swanhubx/swanlab/core/internal/syncer"
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

// 进程退出码约定：0 正常退出（含信号触发的优雅关闭）；1 运行错误；2 用法错误；
// 3 端点已被其它 core 占用。
const (
	exitRunError     = 1
	exitUsageError   = 2
	exitEndpointBusy = 3
)

// 自建端点与收尾参数。
const (
	coreSocketName     = portinfo.CoreSocketFileName
	loopbackAddr       = "127.0.0.1:0"
	shutdownGrace      = 10 * time.Second
	cleanupLockTimeout = 2 * time.Second
)

// endpointProbeTimeout 是端点存活探测的超时上限。
const endpointProbeTimeout = 200 * time.Millisecond

// maxRecvMsgSize 是单条 gRPC 消息的上限。
const maxRecvMsgSize = 4 << 20

// errEndpointBusy 表示端点已被其它 core 占用。
var errEndpointBusy = errors.New("endpoint is already served by a running swanlab-core")

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
	logDir := fs.String("log-dir", "",
		"diagnostic log directory; empty disables file logging")
	detach := fs.Bool("detach", false,
		"detach from parent process; keep running until signal, teardown, or idle timeout")
	idleTimeout := fs.Duration("idle-timeout", 0,
		"idle shutdown timeout when active connections reach zero (e.g. 5m, 2s; 0 disables)")
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
	if *idleTimeout < 0 {
		console.Error("invalid --idle-timeout: must be non-negative")
		return exitUsageError
	}
	if *logDir != "" {
		if err := console.Init(*logDir); err != nil {
			console.Error("failed to open diagnostic log:", err)
			return exitRunError
		}
	}

	var startupLock *portinfo.FileLock
	releaseStartupLock := func() {
		if startupLock != nil {
			_ = startupLock.Release()
			startupLock = nil
		}
	}
	defer releaseStartupLock()

	if *listenAddr == "" && *portFilename == "" {
		if !*detach {
			console.Error("no listen endpoint: pass --listen (manual debug) or --port-filename (SDK startup convention)")
			return exitUsageError
		}
		resolvedPortFile, lock, rtErr := preparePerUserPortFile(context.Background())
		if rtErr != nil {
			console.Error("failed to prepare per-user runtime:", rtErr)
			return exitRunError
		}
		*portFilename = resolvedPortFile
		startupLock = lock
	}

	// port-file 仍指向存活的 core 时不覆盖：覆盖会让原实例失去发现入口。
	if *portFilename != "" {
		if pid, alive := portFileServed(*portFilename); alive {
			console.Errorf("port-file %s is owned by a running swanlab-core (pid %d); refusing to start a second instance", *portFilename, pid)
			return exitEndpointBusy
		}
	}

	// 自建资源记录，退出时持启动锁并核对文件身份后清理自己创建的部分。
	owned := endpointOwnership{pid: os.Getpid()}
	defer func() {
		cleanupEndpoint(&owned, cleanupLockTimeout)
	}()

	ln, err := openEndpoint(*listenAddr, *portFilename)
	if err != nil {
		console.Error("listen failed:", err)
		if errors.Is(err, errEndpointBusy) {
			return exitEndpointBusy
		}
		return exitRunError
	}
	defer func() { _ = ln.Close() }()
	trackListenerOwnership(ln, &owned)

	parentExited, pid, watchErr := watchParentIfAttached(*detach, *parentPID)
	if watchErr != nil {
		console.Error("failed to watch parent process, aborting startup:", watchErr)
		return exitRunError
	}

	// 装配服务核心分层：
	// - server.Controller: 负责 gRPC 进程宿主生命周期与优雅退出控制；
	// - manager.Manager: 负责 Run 领域会话管理与路由注册；
	// - service.CoreService: 负责 gRPC 请求接入与协议映射。
	grpcServer, ctrl, mgr, syncMgr, idleMonitor := assembleServer(*idleTimeout)
	defer func() {
		_ = mgr.Close()
		_ = syncMgr.Close()
	}()

	// listen 与 server 初始化成功后写 port-file，并立即释放启动锁。
	if *portFilename != "" {
		if writeErr := writePortFile(*portFilename, owned.pid, ln.Addr()); writeErr != nil {
			console.Error("failed to write port-file:", writeErr)
			return exitRunError
		}
		owned.portFile = *portFilename
		if pfInfo, statErr := os.Lstat(*portFilename); statErr == nil {
			owned.portFileInfo = pfInfo
		}
	}
	releaseStartupLock()

	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer stop()

	if *detach {
		console.Infof("swanlab-core %s listening on %s (detached, idle timeout %s)", version, ln.Addr(), *idleTimeout)
	} else {
		console.Infof("swanlab-core %s listening on %s (parent pid %d)", version, ln.Addr(), pid)
	}

	serveErr := make(chan error, 1)
	go func() {
		serveErr <- grpcServer.Serve(ln)
	}()
	if idleMonitor != nil {
		idleMonitor.Start()
	}

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
	if serveFailure != nil && !errors.Is(serveFailure, grpc.ErrServerStopped) {
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
		ln, err := listenAutoUnix(sockPath)
		if err == nil {
			return ln, nil
		}
		if errors.Is(err, errEndpointBusy) {
			return nil, err
		}
		console.Warningf("unix listen on %s failed: %v; falling back to loopback tcp %s", sockPath, err, loopbackAddr)
		return net.Listen("tcp", loopbackAddr)
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
		if !strings.HasPrefix(rest, "/") && !strings.HasPrefix(rest, "@") {
			return nil, fmt.Errorf("unix listen address must be absolute or abstract: %q", rest)
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

// listenAutoUnix 在自动模式绑定 UDS。占用者是活的监听则返回 errEndpointBusy；
// 是崩溃残留的 socket 则回收后重试；其它文件不动，返回错误由调用方回退 TCP。
func listenAutoUnix(sockPath string) (net.Listener, error) {
	ln, err := net.Listen("unix", sockPath)
	if err == nil {
		return ln, nil
	}
	if !errors.Is(err, syscall.EADDRINUSE) {
		return nil, err
	}
	if endpointAlive("unix", sockPath) {
		return nil, fmt.Errorf("%w: %s", errEndpointBusy, sockPath)
	}
	info, statErr := os.Lstat(sockPath)
	if statErr != nil || info.Mode()&os.ModeSocket == 0 {
		return nil, err
	}
	if rmErr := os.Remove(sockPath); rmErr != nil {
		return nil, fmt.Errorf("remove stale socket %s: %w", sockPath, rmErr)
	}
	console.Warningf("removed stale core socket %s left by a previous instance", sockPath)
	return net.Listen("unix", sockPath)
}

// endpointAlive 用短超时连接端点，判断是否有进程正在监听。
func endpointAlive(network, addr string) bool {
	conn, err := net.DialTimeout(network, addr, endpointProbeTimeout)
	if err != nil {
		return false
	}
	_ = conn.Close()
	return true
}

// portFileServed 判断 port-file 指向的端点是否仍有 core 在服务，返回其 pid。
// 文件缺失或端点已失效时按可覆盖处理。
func portFileServed(path string) (int, bool) {
	info, err := portinfo.ParseFile(path)
	if err != nil {
		return 0, false
	}
	if info.UnixPath != "" {
		return info.PID, endpointAlive("unix", info.UnixPath)
	}
	return info.PID, endpointAlive("tcp", net.JoinHostPort("127.0.0.1", strconv.Itoa(info.SockPort)))
}

type endpointOwnership struct {
	pid          int
	portFile     string
	portFileInfo os.FileInfo
	socketPath   string
	socketInfo   os.FileInfo
}

func trackListenerOwnership(ln net.Listener, owned *endpointOwnership) {
	if unixLn, ok := ln.(*net.UnixListener); ok {
		// 禁用 Close() 自动 unlink，避免旧实例退出时无锁删掉新实例刚重建的 socket。
		unixLn.SetUnlinkOnClose(false)
	}
	if addr, ok := ln.Addr().(*net.UnixAddr); ok && !strings.HasPrefix(addr.Name, "@") {
		owned.socketPath = addr.Name
		if sockInfo, statErr := os.Lstat(addr.Name); statErr == nil {
			owned.socketInfo = sockInfo
		}
	}
}

// cleanupEndpoint 串行化退出清理与新实例接管：优先在运行时目录持有 core.lock，
// 再分别核对 port-file 与 UDS socket 的文件身份，确认未被新实例替换后才删除。
func cleanupEndpoint(owned *endpointOwnership, lockTimeout time.Duration) {
	if owned == nil || (owned.portFile == "" && owned.socketPath == "") {
		return
	}
	if owned.portFile != "" {
		lockDir := filepath.Dir(owned.portFile)
		if portinfo.ValidateRuntimeDir(lockDir) == nil {
			ctx, cancel := context.WithTimeout(context.Background(), lockTimeout)
			defer cancel()
			lock, err := portinfo.AcquireLock(ctx, lockDir, portinfo.LockOptions{Timeout: lockTimeout})
			if err != nil {
				// 有新实例正在竞争或持锁接管时宁可保留旧文件，交由后续启动流程清理。
				return
			}
			defer func() { _ = lock.Release() }()
		}
		// Linux tmpfs/ext4 在 unlink+新建后可能立刻复用同一 inode；若 port-file 已被新 PID 接管，
		// 必须连同 core.sock 一起保留，不能仅凭 inode 核对误删。
		if runtimeOwnedByOther(owned.portFile, owned.pid) {
			return
		}
	}
	cleanupPortFileChecked(owned.portFile, owned.pid, owned.portFileInfo)
	cleanupSocketChecked(owned.socketPath, owned.socketInfo)
}

func runtimeOwnedByOther(portFile string, selfPID int) bool {
	if portFile == "" || selfPID <= 0 {
		return false
	}
	info, err := portinfo.ParseFile(portFile)
	if err != nil {
		return false
	}
	return info.PID != selfPID
}

func sameFileAtPath(path string, expected os.FileInfo) bool {
	if path == "" || expected == nil {
		return false
	}
	latest, err := os.Lstat(path)
	return err == nil && os.SameFile(expected, latest)
}

func cleanupSocketChecked(path string, expected os.FileInfo) {
	if strings.HasPrefix(path, "@") || !sameFileAtPath(path, expected) {
		return
	}
	// 即使未走 port-file 且 inode 恰好被复用，也不得删除仍有活进程监听的 socket。
	if endpointAlive("unix", path) {
		return
	}
	_ = os.Remove(path)
}

func cleanupPortFileChecked(path string, pid int, expected os.FileInfo) {
	if pid <= 0 || !sameFileAtPath(path, expected) {
		return
	}
	info, err := portinfo.ParseFile(path)
	if err != nil || info.PID != pid {
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

// preparePerUserPortFile 解析并初始化固定 per-user 运行时目录，抢占启动锁后返回 port-file 路径。
func preparePerUserPortFile(ctx context.Context) (string, *portinfo.FileLock, error) {
	runtimeDir, err := portinfo.DefaultRuntimeDir()
	if err != nil {
		return "", nil, err
	}
	return preparePortFileInDir(ctx, runtimeDir)
}

func preparePortFileInDir(ctx context.Context, runtimeDir string) (string, *portinfo.FileLock, error) {
	if err := portinfo.EnsureRuntimeDir(runtimeDir); err != nil {
		return "", nil, err
	}
	if err := portinfo.ValidateSocketPath(runtimeDir); err != nil {
		return "", nil, err
	}
	lock, err := portinfo.AcquireLock(ctx, runtimeDir, portinfo.LockOptions{})
	if err != nil {
		return "", nil, err
	}
	return filepath.Join(runtimeDir, portinfo.CorePortFileName), lock, nil
}

func assembleServer(
	idleTimeout time.Duration,
) (*grpc.Server, *server.Controller, *manager.Manager, *syncer.Manager, *server.IdleMonitor) {
	httpOpts := client.Options{Version: version}
	mgr := manager.New(manager.Options{HTTP: httpOpts})
	syncMgr := syncer.NewManager(syncer.Options{HTTP: httpOpts})

	tracker := server.NewConnTracker()
	grpcServer := grpc.NewServer(
		grpc.MaxRecvMsgSize(maxRecvMsgSize),
		grpc.StatsHandler(tracker),
		grpc.UnaryInterceptor(newPanicRecoveryInterceptor(mgr.Invalidate)),
	)
	ctrl := server.NewController(grpcServer, shutdownGrace)
	var idleMonitor *server.IdleMonitor
	if idleTimeout > 0 {
		idleMonitor = server.NewIdleMonitor(idleTimeout, func() {
			ctrl.Shutdown("idle-timeout")
		})
		tracker.SetMonitor(idleMonitor)
		ctrl.SetIdleMonitor(idleMonitor)
	}
	svc := service.NewCoreService(ctrl, mgr)
	syncSvc := service.NewSyncService(ctrl, syncMgr)
	tracker.SetOnDisconnect(func(connID uint64) {
		svc.HandleConnDisconnect(connID)
		syncSvc.HandleConnDisconnect(connID)
	})
	// 关闭顺序：并发等待断连清理与活跃会话排空（共享 shutdownGrace），随后停止 gRPC。
	ctrl.SetOnShutdown(func(ctx context.Context) error {
		var (
			waitRunErr      error
			waitSyncErr     error
			shutdownRunErr  error
			shutdownSyncErr error
			wg              sync.WaitGroup
		)
		wg.Add(3)
		go func() {
			defer wg.Done()
			waitRunErr = svc.WaitDisconnectCleanups(ctx)
		}()
		go func() {
			defer wg.Done()
			waitSyncErr = syncSvc.WaitDisconnectCleanups(ctx)
		}()
		go func() {
			defer wg.Done()
			shutdownSyncErr = syncMgr.Shutdown(ctx)
		}()
		shutdownRunErr = mgr.Shutdown(ctx)
		wg.Wait()
		return errors.Join(waitRunErr, waitSyncErr, shutdownRunErr, shutdownSyncErr)
	})
	svc.Register(grpcServer)
	syncSvc.Register(grpcServer)
	return grpcServer, ctrl, mgr, syncMgr, idleMonitor
}

type sessionRequest interface {
	GetSessionId() string
}

// newPanicRecoveryInterceptor 防止单个 RPC handler 的 panic 终止整个 core 进程，
// 在恢复后立即使异常会话失效，并按方法名对完整调用栈输出做去重节流。
func newPanicRecoveryInterceptor(invalidate func(string)) grpc.UnaryServerInterceptor {
	var loggedStack sync.Map
	return func(ctx context.Context, req any, info *grpc.UnaryServerInfo, handler grpc.UnaryHandler) (resp any, err error) {
		defer func() {
			if recovered := recover(); recovered != nil {
				method := "<unknown>"
				if info != nil && info.FullMethod != "" {
					method = info.FullMethod
				}
				if _, seen := loggedStack.LoadOrStore(method, struct{}{}); !seen {
					console.Errorf("panic in gRPC %s: %v\n%s", method, recovered, debug.Stack())
				} else {
					console.Errorf("panic in gRPC %s: %v (stack suppressed)", method, recovered)
				}
				if invalidate != nil {
					if sr, ok := req.(sessionRequest); ok {
						if sid := sr.GetSessionId(); sid != "" {
							invalidate(sid)
						}
					}
				}
				resp = nil
				err = status.Errorf(codes.Internal, "internal panic: %v", recovered)
			}
		}()
		return handler(ctx, req)
	}
}

func watchParentIfAttached(detach bool, parentPID int) (<-chan struct{}, int, error) {
	if detach {
		return nil, 0, nil
	}
	pid := parentPID
	if pid <= 0 {
		pid = os.Getppid()
	}
	ch, err := process.NotifyOnParentExit(pid)
	if err != nil {
		return nil, 0, err
	}
	return ch, pid, nil
}

func writePortFile(path string, pid int, addr net.Addr) error {
	info := portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: pid}
	switch typed := addr.(type) {
	case *net.UnixAddr:
		info.UnixPath = typed.Name
	case *net.TCPAddr:
		info.SockPort = typed.Port
	default:
		return fmt.Errorf("unrecognized listener address type: %v", addr)
	}
	return portinfo.WriteFile(path, &info)
}

// envInt 解析整型环境变量，缺失或非法时返回 0。
func envInt(name string) int {
	v, err := strconv.Atoi(strings.TrimSpace(os.Getenv(name)))
	if err != nil {
		return 0
	}
	return v
}
