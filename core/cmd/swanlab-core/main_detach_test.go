package main

import (
	"context"
	"errors"
	"fmt"
	"io"
	"net"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/pkg/portinfo"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	settingsv1 "github.com/swanhubx/swanlab/core/proto/swanlab/settings/core/v1"
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
	// 即使无 port-file（如 --listen 调试场景）且系统复用了同一 inode，也不得删除活的监听 socket。
	if runtime.GOOS != "windows" {
		cleanupEndpoint(&endpointOwnership{
			pid:        1001,
			socketPath: oldSockPath,
			socketInfo: oldSockInfo,
		}, time.Second)
		if !endpointAlive("unix", oldSockPath) {
			t.Fatal("live unix socket was removed without port-file")
		}
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

func TestRunDetachedDisconnectCleansSessionAndExitsIdle(t *testing.T) {
	dir := shortRuntimeDir(t)
	portFile := filepath.Join(dir, portinfo.CorePortFileName)
	runDir := t.TempDir()
	const runID = "detach_disconnect_run"

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

	startResp, startErr := client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: &settingsv1.CoreSettings{
			RunId:          runID,
			RunDir:         runDir,
			Mode:           settingsv1.CoreMode_CORE_MODE_LOCAL,
			RecordBatch:    16,
			RecordInterval: 0.05,
			SaveSize:       1024,
			SaveSplit:      1024,
			SavePart:       1024,
			SaveBatch:      16,
		},
		StartRecord: &runv1.StartRecord{Id: runID, Project: "smoke"},
	})
	if startErr != nil || !startResp.GetSuccess() {
		_ = conn.Close()
		t.Fatalf("DeliverRunStart resp=%v err=%v", startResp, startErr)
	}

	if _, upErr := client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{
		SessionId: startResp.GetSessionId(),
		Data: []*datav1.ScalarRecord{{
			Key:   "loss",
			Step:  1,
			Type:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
			Value: &datav1.ScalarValue{Number: 0.25},
		}},
	}); upErr != nil {
		_ = conn.Close()
		t.Fatalf("UpsertScalars: %v", upErr)
	}

	// 不调 finish/confirm 直接断开连接：Core 须先排空落盘该会话，再在空闲到期后退出。
	if closeErr := conn.Close(); closeErr != nil {
		t.Fatalf("conn.Close: %v", closeErr)
	}

	select {
	case code := <-doneCh:
		if code != 0 {
			t.Fatalf("core exit code = %d, want 0", code)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("detached core did not exit after disconnect and idle timeout")
	}

	runFile := filepath.Join(runDir, "run-"+runID+".swanlab")
	reader, openErr := store.Open(runFile)
	if openErr != nil {
		t.Fatalf("store.Open(%q): %v", runFile, openErr)
	}
	defer func() { _ = reader.Close() }()

	var hasScalar, hasFinish bool
	for {
		payload, scanErr := reader.Scan()
		if errors.Is(scanErr, io.EOF) {
			break
		}
		if scanErr != nil {
			t.Fatalf("reader.Scan: %v", scanErr)
		}
		rec := &recordv1.Record{}
		if unmarshalErr := proto.Unmarshal(payload, rec); unmarshalErr != nil {
			t.Fatalf("proto.Unmarshal: %v", unmarshalErr)
		}
		if rec.GetScalar() != nil && rec.GetScalar().GetKey() == "loss" {
			hasScalar = true
		}
		if rec.GetNum() == record.NumFinish {
			hasFinish = true
		}
	}
	if !hasScalar {
		t.Fatal("disconnect cleanup must flush pending scalar record to local store")
	}
	if hasFinish {
		t.Fatal("disconnect cleanup must not synthesize a FinishRecord")
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

func TestRunDegradesWhenDiagnosticLogOpenFails(t *testing.T) {
	dir := shortRuntimeDir(t)
	portFile := filepath.Join(dir, portinfo.CorePortFileName)

	// 将 debug-core.log 占成目录，模拟日志路径不可写/不可打开的真实失败场景。
	badLogDir := t.TempDir()
	if err := os.Mkdir(filepath.Join(badLogDir, "debug-core.log"), 0o700); err != nil {
		t.Fatalf("mkdir debug-core.log: %v", err)
	}

	r, w, err := os.Pipe()
	if err != nil {
		t.Fatalf("os.Pipe: %v", err)
	}
	oldStderr := os.Stderr
	os.Stderr = w
	stderrCh := make(chan string, 1)
	go func() {
		b, _ := io.ReadAll(r)
		stderrCh <- string(b)
	}()
	restored := false
	restoreStderr := func() string {
		if !restored {
			restored = true
			os.Stderr = oldStderr
			_ = w.Close()
		}
		out := <-stderrCh
		_ = r.Close()
		return out
	}
	defer func() {
		if !restored {
			os.Stderr = oldStderr
			_ = w.Close()
			<-stderrCh
			_ = r.Close()
		}
	}()

	doneCh := make(chan int, 1)
	go func() {
		doneCh <- run([]string{"--port-filename", portFile, "--detach", "--idle-timeout", "5s", "--log-dir", badLogDir})
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
	if _, tdErr := client.TeardownService(ctx, &corev1.TeardownServiceRequest{}); tdErr != nil {
		t.Fatalf("TeardownService: %v", tdErr)
	}

	select {
	case code := <-doneCh:
		if code != 0 {
			t.Fatalf("core exit code = %d, want 0", code)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("core did not exit after TeardownService")
	}

	stderrOut := restoreStderr()
	if !strings.Contains(stderrOut, "failed to open diagnostic log") {
		t.Fatalf("expected stderr warning for diagnostic log failure, got %q", stderrOut)
	}
}
