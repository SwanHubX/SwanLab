// Command swanlab-core 是 SwanLab Go core 的进程入口。
//
// 当前为发布脚手架阶段：仅建立进程生命周期骨架——版本上报（--version）、
// 端点监听、父进程退出监控与信号处理；gRPC 服务端（CoreService /
// CoreSyncService / ProbeService）在后续迭代中接入，届时替换 serve 循环。
//
// 端点约定：
//
//	--listen unix:///path/to/uds   Linux/macOS 进程内通信（默认路径由 Python SDK 分配）
//	--listen tcp://127.0.0.1:port  Windows 回环地址（uds 不可用，named pipe 支持后续提供）
//
// 生命周期：父进程退出（process 包监控）或收到 SIGINT/SIGTERM 时优雅退出，
// 防止 Python SDK 崩溃后 core 沦为孤儿进程。
package main

import (
	"context"
	"errors"
	"flag"
	"fmt"
	"net"
	"os"
	"os/signal"
	"runtime"
	"strconv"
	"strings"
	"syscall"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/process"
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

func main() {
	os.Exit(run(os.Args[1:]))
}

func run(args []string) int {
	fs := flag.NewFlagSet("swanlab-core", flag.ContinueOnError)
	printVersion := fs.Bool("version", false, "打印版本信息后退出")
	listenAddr := fs.String("listen", os.Getenv(envListenAddr),
		"监听端点，格式 unix://<uds路径> 或 tcp://<地址:端口>；Windows 仅支持 tcp:// 回环地址")
	parentPID := fs.Int("parent-pid", envInt(envParentPID),
		"预期父进程 PID，父进程退出时 core 随之退出；未指定时取启动瞬间的实际父进程")
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
	if *listenAddr == "" {
		console.Error("未指定监听端点：通过 --listen 或环境变量 " + envListenAddr + " 传入")
		return exitUsageError
	}

	ln, err := listen(*listenAddr)
	if err != nil {
		console.Error("监听失败:", err)
		return exitRunError
	}
	defer ln.Close()

	// 父进程监控：显式传入的 PID 优先（Python SDK 启动约定），未传时回退为
	// 监控启动瞬间的实际父进程（本地终端运行场景）。监控建立失败按约定终止启动。
	pid := *parentPID
	if pid <= 0 {
		pid = os.Getppid()
	}
	parentExited, err := process.NotifyOnParentExit(pid)
	if err != nil {
		console.Error("父进程监控建立失败，终止启动:", err)
		return exitRunError
	}

	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer stop()

	console.Infof("swanlab-core %s listening on %s (parent pid %d)", version, *listenAddr, pid)

	serveErr := make(chan error, 1)
	go func() {
		serveErr <- serve(ln)
	}()

	select {
	case <-ctx.Done():
		console.Info("收到退出信号，正在关闭")
	case <-parentExited:
		console.Warning("父进程已退出，core 随之退出")
	case err := <-serveErr:
		if err != nil {
			console.Error("监听异常退出:", err)
			return exitRunError
		}
	}
	return 0
}

// listen 按协议前缀创建监听器。uds 仅在非 Windows 平台可用；Windows 使用
// TCP 回环地址兜底（named pipe 接入后在此分支扩展）。
func listen(addr string) (net.Listener, error) {
	scheme, rest, ok := strings.Cut(addr, "://")
	if !ok {
		return nil, fmt.Errorf("监听端点缺少协议前缀（unix:// 或 tcp://）: %s", addr)
	}
	switch scheme {
	case "unix":
		if runtime.GOOS == "windows" {
			return nil, errors.New("windows 平台不支持 unix:// 端点，请使用 tcp://127.0.0.1:<端口>")
		}
		return net.Listen("unix", rest)
	case "tcp":
		return net.Listen("tcp", rest)
	default:
		return nil, fmt.Errorf("不支持的监听协议 %q（仅 unix:// 或 tcp://）", scheme)
	}
}

// serve 接受连接后立即关闭。脚手架阶段仅验证端点连通性，
// gRPC 服务端就绪后由此接入 serve 逻辑。
func serve(ln net.Listener) error {
	for {
		conn, err := ln.Accept()
		if err != nil {
			if errors.Is(err, net.ErrClosed) {
				return nil
			}
			return err
		}
		_ = conn.Close()
	}
}

// envInt 解析整型环境变量，缺失或非法时返回 0。
func envInt(name string) int {
	v, err := strconv.Atoi(strings.TrimSpace(os.Getenv(name)))
	if err != nil {
		return 0
	}
	return v
}
