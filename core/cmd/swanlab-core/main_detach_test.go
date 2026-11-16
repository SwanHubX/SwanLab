package main

import (
	"context"
	"errors"
	"fmt"
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
