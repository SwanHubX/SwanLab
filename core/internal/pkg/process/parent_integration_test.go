package process

import (
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"syscall"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// parentExitHelperMode 仅用于让同一个测试二进制分别扮演辅助父进程和被监控子进程。
// 测试进程 -> parent helper -> monitored child 的三级结构可以真实验证父进程退出行为。
const parentExitHelperMode = "SWANLAB_PARENT_EXIT_HELPER"

func TestNotifyOnParentExitNotifiesCaller(t *testing.T) {
	switch os.Getenv(parentExitHelperMode) {
	case "parent":
		runParentHelper()
	case "child":
		runMonitoredChildHelper()
	}

	helperDir := t.TempDir()
	cmd := exec.Command(os.Args[0], "-test.run=^TestNotifyOnParentExitNotifiesCaller$")
	cmd.Env = append(os.Environ(), parentExitHelperMode+"=parent", "SWANLAB_PARENT_EXIT_DIR="+helperDir)
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		t.Fatal(err)
	}
	cmd.Stderr = os.Stderr
	if err = cmd.Start(); err != nil {
		t.Fatal(err)
	}

	closed := make(chan error, 1)
	go func() {
		_, copyErr := io.Copy(io.Discard, stdout)
		closed <- copyErr
	}()

	select {
	case waitErr := <-closed:
		if waitErr != nil {
			t.Fatalf("wait for monitored process: %v", waitErr)
		}
	case <-time.After(10 * time.Second):
		killHelperProcess(t, filepath.Join(helperDir, "child.pid"))
		t.Fatal("caller did not exit after receiving the parent exit notification")
	}

	if err = cmd.Wait(); err != nil {
		t.Fatalf("parent helper failed: %v", err)
	}
	data, err := os.ReadFile(filepath.Join(helperDir, "debug-core.log"))
	if err != nil {
		t.Fatalf("read console log: %v", err)
	}
	if !strings.Contains(string(data), "Monitoring parent process") {
		t.Fatalf("console log does not contain monitor setup: %q", data)
	}
}

func runParentHelper() {
	helperDir := os.Getenv("SWANLAB_PARENT_EXIT_DIR")
	cmd := exec.Command(os.Args[0], "-test.run=^TestNotifyOnParentExitNotifiesCaller$")
	cmd.Env = append(
		os.Environ(),
		parentExitHelperMode+"=child",
		"SWANLAB_PARENT_EXIT_DIR="+helperDir,
		"SWANLAB_DEBUG=true",
	)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Start(); err != nil {
		fmt.Fprintf(os.Stderr, "start monitored child: %v\n", err)
		os.Exit(2)
	}

	readyPath := filepath.Join(helperDir, "ready")
	deadline := time.Now().Add(5 * time.Second)
	for time.Now().Before(deadline) {
		if _, err := os.Stat(readyPath); err == nil {
			os.Exit(0)
		}
		time.Sleep(10 * time.Millisecond)
	}
	fmt.Fprintln(os.Stderr, "monitored child did not become ready")
	os.Exit(2)
}

func runMonitoredChildHelper() {
	helperDir := os.Getenv("SWANLAB_PARENT_EXIT_DIR")
	if err := console.Init(helperDir); err != nil {
		fmt.Fprintf(os.Stderr, "initialize console: %v\n", err)
		os.Exit(2)
	}
	parentExited, err := NotifyOnParentExit(os.Getppid())
	if err != nil {
		fmt.Fprintf(os.Stderr, "monitor parent: %v\n", err)
		os.Exit(2)
	}
	if err := os.WriteFile(filepath.Join(helperDir, "child.pid"), []byte(strconv.Itoa(os.Getpid())), 0o600); err != nil {
		fmt.Fprintf(os.Stderr, "write child PID: %v\n", err)
		os.Exit(2)
	}
	if err := os.WriteFile(filepath.Join(helperDir, "ready"), nil, 0o600); err != nil {
		fmt.Fprintf(os.Stderr, "mark child ready: %v\n", err)
		os.Exit(2)
	}

	// process 包只负责发出通知；由调用方决定收到通知后的退出方式。
	<-parentExited
	os.Exit(0)
}

// TestNotifyOnParentExitSurvivesSpawnThreadExit 验证进程级父死兜底（G5）：
//
// 场景：父进程从非主线程 spawn 被监控子进程，该线程先退出（触发一次 PDEATHSIG
// 伪通知，必须被 PPID 复查过滤），父进程本体随后才退出。PDEATHSIG 绑定的是
// spawn 线程，伪通知之后不再有信号，子进程只能依赖 pidfd/PPID 轮询兜底感知
// 父进程退出，且不能因伪通知提前退出。
//
// 三级结构与 TestNotifyOnParentExitNotifiesCaller 一致：
// 测试进程 -> parent helper（从 LockOSThread 的 goroutine spawn child）-> monitored child。
func TestNotifyOnParentExitSurvivesSpawnThreadExit(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PDEATHSIG spawn-thread semantics are Linux-specific")
	}

	switch os.Getenv(parentExitHelperMode) {
	case "thread-parent":
		runThreadParentHelper()
	case "thread-child":
		runThreadMonitoredChildHelper()
	}

	helperDir := t.TempDir()
	cmd := exec.Command(os.Args[0], "-test.run=^TestNotifyOnParentExitSurvivesSpawnThreadExit$")
	cmd.Env = append(
		os.Environ(),
		parentExitHelperMode+"=thread-parent",
		"SWANLAB_PARENT_EXIT_DIR="+helperDir,
	)
	cmd.Stderr = os.Stderr
	if err := cmd.Start(); err != nil {
		t.Fatal(err)
	}

	// 1. child 就绪且 spawn 线程已退出后，child 不得因 PDEATHSIG 伪通知提前退出。
	if err := waitForMarker(filepath.Join(helperDir, "thread-dead"), 5*time.Second); err != nil {
		killHelperProcess(t, filepath.Join(helperDir, "child.pid"))
		t.Fatalf("wait for spawn thread exit: %v", err)
	}
	childPID, err := readPIDFile(filepath.Join(helperDir, "child.pid"))
	if err != nil {
		t.Fatalf("read child pid: %v", err)
	}
	if !processAlive(childPID) {
		killHelperProcess(t, filepath.Join(helperDir, "child.pid"))
		t.Fatal("child exited prematurely after spawn thread death (spurious PDEATHSIG not filtered)")
	}

	// 2. parent helper 自行退出（thread-parent 模式在 thread-dead 后 os.Exit(0)），
	//    child 必须在兜底监控预算内感知并退出。
	if err := cmd.Wait(); err != nil {
		t.Fatalf("thread parent helper failed: %v", err)
	}
	// parent 已退出但可能尚未被内核重新托管，轮询等待 child 退出。
	deadline := time.Now().Add(10 * time.Second)
	for time.Now().Before(deadline) {
		if !processAlive(childPID) {
			return // child 正确随父进程退出
		}
		time.Sleep(50 * time.Millisecond)
	}
	t.Fatal("child did not exit after parent process death (process-level fallback failed)")
}

// runThreadParentHelper 在专用 OS 线程上 spawn child，随后结束该 goroutine
// 使线程销毁（等价于"Python 从非主线程 spawn 且该线程先退出"），确认 child
// 存活后退出自身，验证 child 的进程级父死监控。
func runThreadParentHelper() {
	helperDir := os.Getenv("SWANLAB_PARENT_EXIT_DIR")
	type result struct {
		pid int
		err error
	}
	ch := make(chan result, 1)
	go func() {
		// 锁定 OS 线程：本 goroutine 返回时线程被 runtime 终止，
		// 内核随之向 child 发送一次 PDEATHSIG 伪通知。
		runtime.LockOSThread()
		cmd := exec.Command(os.Args[0], "-test.run=^TestNotifyOnParentExitSurvivesSpawnThreadExit$")
		cmd.Env = append(
			os.Environ(),
			parentExitHelperMode+"=thread-child",
			"SWANLAB_PARENT_EXIT_DIR="+helperDir,
		)
		cmd.Stderr = os.Stderr
		if err := cmd.Start(); err != nil {
			ch <- result{err: err}
			return
		}
		ch <- result{pid: cmd.Process.Pid}
		// goroutine 返回 -> 锁定线程销毁 -> PDEATHSIG 伪通知
	}()
	res := <-ch
	if res.err != nil {
		fmt.Fprintf(os.Stderr, "start monitored child: %v\n", res.err)
		os.Exit(2)
	}
	if err := os.WriteFile(filepath.Join(helperDir, "child.pid"), []byte(strconv.Itoa(res.pid)), 0o600); err != nil {
		fmt.Fprintf(os.Stderr, "write child PID: %v\n", err)
		os.Exit(2)
	}
	// 等 child 建立监控（ready 标记），再给线程销毁留出传播时间。
	if err := waitForMarkerFile(filepath.Join(helperDir, "ready"), 5*time.Second); err != nil {
		fmt.Fprintf(os.Stderr, "monitored child did not become ready: %v\n", err)
		os.Exit(2)
	}
	time.Sleep(300 * time.Millisecond)
	if err := os.WriteFile(filepath.Join(helperDir, "thread-dead"), nil, 0o600); err != nil {
		fmt.Fprintf(os.Stderr, "mark thread dead: %v\n", err)
		os.Exit(2)
	}
	// 覆盖至少一个轮询周期，证明 child 在线程死亡 + 伪通知后仍然存活。
	time.Sleep(1500 * time.Millisecond)
	os.Exit(0)
}

// runThreadMonitoredChildHelper 与 runMonitoredChildHelper 类似，
// 但在父进程退出后写 exited 标记再退出，供测试区分"随父退出"与"提前退出"。
func runThreadMonitoredChildHelper() {
	helperDir := os.Getenv("SWANLAB_PARENT_EXIT_DIR")
	parentExited, err := NotifyOnParentExit(os.Getppid())
	if err != nil {
		fmt.Fprintf(os.Stderr, "monitor parent: %v\n", err)
		os.Exit(2)
	}
	if err := os.WriteFile(filepath.Join(helperDir, "ready"), nil, 0o600); err != nil {
		fmt.Fprintf(os.Stderr, "mark child ready: %v\n", err)
		os.Exit(2)
	}
	// process 包只负责发出通知；由调用方决定收到通知后的退出方式。
	<-parentExited
	_ = os.WriteFile(filepath.Join(helperDir, "exited"), nil, 0o600)
	os.Exit(0)
}

func readPIDFile(path string) (int, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return 0, err
	}
	return strconv.Atoi(string(data))
}

// processAlive 通过 signal 0 探测进程是否存在（不含僵尸态判断，
// 测试场景中进程由测试进程树的子孙构成，不会被长期悬挂为僵尸）。
func processAlive(pid int) bool {
	proc, err := os.FindProcess(pid)
	if err != nil {
		return false
	}
	return proc.Signal(syscall.Signal(0)) == nil
}

func waitForMarker(path string, timeout time.Duration) error {
	return waitForMarkerFile(path, timeout)
}

func waitForMarkerFile(path string, timeout time.Duration) error {
	deadline := time.Now().Add(timeout)
	for time.Now().Before(deadline) {
		if _, err := os.Stat(path); err == nil {
			return nil
		}
		time.Sleep(10 * time.Millisecond)
	}
	return fmt.Errorf("marker %s not found within %s", path, timeout)
}

func killHelperProcess(t *testing.T, pidPath string) {
	data, err := os.ReadFile(pidPath)
	if err != nil {
		t.Logf("read helper PID for cleanup: %v", err)
		return
	}
	pid, err := strconv.Atoi(string(data))
	if err != nil {
		t.Logf("parse helper PID for cleanup: %v", err)
		return
	}
	process, err := os.FindProcess(pid)
	if err == nil {
		_ = process.Kill()
	}
}
