//go:build windows

package process

import (
	"fmt"

	"golang.org/x/sys/windows"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

func notifyOnParentExit(parentPID int) (<-chan struct{}, error) {
	// Windows 中父进程退出后 PPID 不会像 Unix reparent 那样改变；持有带 SYNCHRONIZE
	// 权限的进程句柄并等待它变为 signaled，才能可靠判断原父进程已经结束。
	parent, err := windows.OpenProcess(windows.SYNCHRONIZE, false, uint32(parentPID))
	if err != nil {
		return nil, fmt.Errorf("open parent process %d: %w", parentPID, err)
	}
	if err := checkParent(parentPID); err != nil {
		_ = windows.CloseHandle(parent)
		return nil, err
	}

	parentExited := make(chan struct{})
	go func() {
		// 句柄由等待 goroutine 独占，直至父进程退出或等待失败后关闭。
		defer windows.CloseHandle(parent)
		status, err := windows.WaitForSingleObject(parent, windows.INFINITE)
		if err != nil {
			console.Error("Failed to wait for parent process:", err)
		} else if status != windows.WAIT_OBJECT_0 {
			console.Errorf("Unexpected parent process wait status: %d", status)
		} else {
			console.Debugf("Parent process %d exited", parentPID)
		}
		// 等待失败时同样通知调用方，避免监控失效后留下孤儿进程。
		close(parentExited)
	}()
	return parentExited, nil
}
