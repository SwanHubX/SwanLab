package process

import (
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
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
