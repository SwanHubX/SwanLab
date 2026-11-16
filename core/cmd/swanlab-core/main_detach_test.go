package main

import (
	"context"
	"errors"
	"fmt"
	"net"
	"os"
	"path/filepath"
	"runtime"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/credentials/insecure"

	"github.com/swanhubx/swanlab/core/internal/pkg/portinfo"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
)

func shortRuntimeDir(t *testing.T) string {
	t.Helper()
	if runtime.GOOS == "windows" {
		return t.TempDir()
	}
	dir, err := os.MkdirTemp("/tmp", "slc")
	if err != nil {
		t.Fatalf("mkdtemp: %v", err)
	}
	t.Cleanup(func() { _ = os.RemoveAll(dir) })
	return dir
}

func waitForPortInfo(t *testing.T, portFile string, doneCh <-chan int) portinfo.Info {
	t.Helper()
	deadline := time.NewTimer(5 * time.Second)
	defer deadline.Stop()
	ticker := time.NewTicker(10 * time.Millisecond)
	defer ticker.Stop()

	for {
		info, err := portinfo.ParseFile(portFile)
		if err == nil {
			return info
		}
		select {
		case code := <-doneCh:
			t.Fatalf("core exited early with code %d before writing %s", code, portFile)
		case <-deadline.C:
			t.Fatalf("timed out waiting for port-file %s", portFile)
		case <-ticker.C:
		}
	}
}

func grpcTarget(info portinfo.Info) string {
	if info.UnixPath != "" {
		return "unix://" + info.UnixPath
	}
	return fmt.Sprintf("127.0.0.1:%d", info.SockPort)
}

func TestRunRejectsNegativeIdleTimeout(t *testing.T) {
	code := run([]string{"--listen", "tcp://127.0.0.1:0", "--idle-timeout", "-1s"})
	if code != exitUsageError {
		t.Fatalf("exit code = %d, want %d", code, exitUsageError)
	}
}

func TestRunRejectsMissingEndpointWithoutDetach(t *testing.T) {
	code := run(nil)
	if code != exitUsageError {
		t.Fatalf("exit code = %d, want %d", code, exitUsageError)
	}
}

func TestPreparePortFileInDirAcquiresAndReleasesLock(t *testing.T) {
	dir := filepath.Join(shortRuntimeDir(t), portinfo.RuntimeDirPrefix+"unit")
	portFile, lock, err := preparePortFileInDir(context.Background(), dir)
	if err != nil {
		t.Fatalf("preparePortFileInDir: %v", err)
	}
	if portFile != filepath.Join(dir, portinfo.CorePortFileName) {
		t.Fatalf("portFile = %q, want %q", portFile, filepath.Join(dir, portinfo.CorePortFileName))
	}
	if lock == nil {
		t.Fatal("expected non-nil startup lock")
	}
	if relErr := lock.Release(); relErr != nil {
		t.Fatalf("Release: %v", relErr)
	}
}

func TestRunDetachedIdleTimeoutAfterClientDisconnect(t *testing.T) {
	dir := shortRuntimeDir(t)
	portFile := filepath.Join(dir, portinfo.CorePortFileName)

	doneCh := make(chan int, 1)
	go func() {
		doneCh <- run([]string{"--port-filename", portFile, "--detach", "--idle-timeout", "100ms"})
	}()

	info := waitForPortInfo(t, portFile, doneCh)
	conn, err := grpc.NewClient(grpcTarget(info), grpc.WithTransportCredentials(insecure.NewCredentials()))
	if err != nil {
		t.Fatalf("grpc.NewClient: %v", err)
	}

	client := corev1.NewCoreServiceClient(conn)
	ctx, cancel := context.WithTimeout(context.Background(), 3*time.Second)
	defer cancel()
	if _, spinErr := client.SpinupService(ctx, &corev1.SpinupServiceRequest{}); spinErr != nil {
		_ = conn.Close()
		t.Fatalf("SpinupService: %v", spinErr)
	}

	// 断开全部连接后，空闲计时器应在超时后触发优雅退出。
	if closeErr := conn.Close(); closeErr != nil {
		t.Fatalf("conn.Close: %v", closeErr)
	}

	select {
	case code := <-doneCh:
		if code != 0 {
			t.Fatalf("core exit code = %d, want 0", code)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("detached core did not exit after idle timeout")
	}

	if _, statErr := os.Stat(portFile); !errors.Is(statErr, os.ErrNotExist) {
		t.Fatalf("port-file still exists after idle exit: %v", statErr)
	}
	if info.UnixPath != "" {
		if _, statErr := os.Stat(info.UnixPath); !errors.Is(statErr, os.ErrNotExist) {
			t.Fatalf("socket file still exists after idle exit: %v", statErr)
		}
	}
}

func TestRunDetachedExitsBusyWhenPortFileServed(t *testing.T) {
	dir := shortRuntimeDir(t)
	portFile := filepath.Join(dir, portinfo.CorePortFileName)

	doneCh := make(chan int, 1)
	go func() {
		doneCh <- run([]string{"--port-filename", portFile, "--detach", "--idle-timeout", "5s"})
	}()

	info := waitForPortInfo(t, portFile, doneCh)
	conn, err := grpc.NewClient(grpcTarget(info), grpc.WithTransportCredentials(insecure.NewCredentials()))
	if err != nil {
		t.Fatalf("grpc.NewClient: %v", err)
	}
	defer func() { _ = conn.Close() }()

	client := corev1.NewCoreServiceClient(conn)
	ctx, cancel := context.WithTimeout(context.Background(), 3*time.Second)
	defer cancel()
	if _, spinErr := client.SpinupService(ctx, &corev1.SpinupServiceRequest{}); spinErr != nil {
		t.Fatalf("SpinupService: %v", spinErr)
	}

	secondCode := run([]string{"--port-filename", portFile, "--detach", "--idle-timeout", "5s"})
	if secondCode != exitEndpointBusy {
		t.Fatalf("second core exit code = %d, want %d", secondCode, exitEndpointBusy)
	}

	if _, tdErr := client.TeardownService(ctx, &corev1.TeardownServiceRequest{}); tdErr != nil {
		t.Fatalf("TeardownService: %v", tdErr)
	}
	select {
	case code := <-doneCh:
		if code != 0 {
			t.Fatalf("first core exit code = %d, want 0", code)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("first core did not exit after TeardownService")
	}
}

func TestCleanupEndpointSkipsReplacedFiles(t *testing.T) {
	dir := shortRuntimeDir(t)
	if ensureErr := portinfo.EnsureRuntimeDir(dir); ensureErr != nil {
		t.Fatalf("EnsureRuntimeDir: %v", ensureErr)
	}
	portFile := filepath.Join(dir, portinfo.CorePortFileName)

	oldPort := portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: 1001, SockPort: 12345}
	var oldSockPath string
	var oldSockInfo os.FileInfo
	var newLn net.Listener

	if runtime.GOOS != "windows" {
		oldSockPath = filepath.Join(dir, coreSocketName)
		oldLn, listenErr := net.Listen("unix", oldSockPath)
		if listenErr != nil {
			t.Fatalf("old unix listen: %v", listenErr)
		}
		oldLn.(*net.UnixListener).SetUnlinkOnClose(false)
		var statErr error
		oldSockInfo, statErr = os.Lstat(oldSockPath)
		if statErr != nil {
			t.Fatal(statErr)
		}
		_ = oldLn.Close()
		// Close() 不得自动 unlink 残留 socket，由统一清理路径负责。
		if _, afterCloseErr := os.Lstat(oldSockPath); afterCloseErr != nil {
			t.Fatalf("SetUnlinkOnClose(false) still removed socket: %v", afterCloseErr)
		}
		oldPort = portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: 1001, UnixPath: oldSockPath}
	}

	if writeErr := portinfo.WriteFile(portFile, &oldPort); writeErr != nil {
		t.Fatal(writeErr)
	}
	oldPortInfo, statErr := os.Lstat(portFile)
	if statErr != nil {
		t.Fatal(statErr)
	}

	// 模拟新实例接管并重建 socket 与 port-file。
	newPort := portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: 2002, SockPort: 23456}
	if runtime.GOOS != "windows" {
		if rmErr := os.Remove(oldSockPath); rmErr != nil {
			t.Fatal(rmErr)
		}
		var listenErr error
		newLn, listenErr = net.Listen("unix", oldSockPath)
		if listenErr != nil {
			t.Fatalf("new unix listen: %v", listenErr)
		}
		defer func() { _ = newLn.Close() }()
		newPort = portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: 2002, UnixPath: oldSockPath}
	}
	if writeErr := portinfo.WriteFile(portFile, &newPort); writeErr != nil {
		t.Fatal(writeErr)
	}

	cleanupEndpoint(&endpointOwnership{
		pid:          1001,
		portFile:     portFile,
		portFileInfo: oldPortInfo,
		socketPath:   oldSockPath,
		socketInfo:   oldSockInfo,
	}, time.Second)

	got, parseErr := portinfo.ParseFile(portFile)
	if parseErr != nil {
		t.Fatalf("new port-file was removed or corrupted: %v", parseErr)
	}
	if got.PID != 2002 {
		t.Fatalf("port-file PID = %d, want 2002", got.PID)
	}
	if runtime.GOOS != "windows" && !endpointAlive("unix", oldSockPath) {
		t.Fatal("new unix socket was removed by old instance cleanup")
	}
}

func TestCleanupEndpointSkipsWhenLockHeld(t *testing.T) {
	dir := shortRuntimeDir(t)
	if ensureErr := portinfo.EnsureRuntimeDir(dir); ensureErr != nil {
		t.Fatalf("EnsureRuntimeDir: %v", ensureErr)
	}
	portFile := filepath.Join(dir, portinfo.CorePortFileName)
	info := portinfo.Info{Protocol: portinfo.ProtocolVersion, PID: 1001, SockPort: 12345}
	if writeErr := portinfo.WriteFile(portFile, &info); writeErr != nil {
		t.Fatal(writeErr)
	}
	portInfo, statErr := os.Lstat(portFile)
	if statErr != nil {
		t.Fatal(statErr)
	}

	lock, acquired, lockErr := portinfo.TryAcquireLock(dir, portinfo.LockOptions{})
	if lockErr != nil || !acquired {
		t.Fatalf("TryAcquireLock: acquired=%v err=%v", acquired, lockErr)
	}
	defer func() { _ = lock.Release() }()

	cleanupEndpoint(&endpointOwnership{
		pid:          1001,
		portFile:     portFile,
		portFileInfo: portInfo,
	}, 20*time.Millisecond)

	if _, pfErr := os.Stat(portFile); pfErr != nil {
		t.Fatalf("port-file must not be deleted while startup lock is held: %v", pfErr)
	}
}
