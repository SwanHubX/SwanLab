package process

import (
	"errors"
	"fmt"
	"os"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// ErrParentExited 表示当前父进程已经不是调用方传入的预期父进程。
// 这通常意味着父进程在监控建立前已经退出，当前进程已被系统重新托管。
var ErrParentExited = errors.New("parent process has exited")

// NotifyOnParentExit 建立父进程退出监控，并在 parentPID 退出时关闭返回的 channel。
//
// parentPID 必须是调用方在启动阶段捕获的预期父进程 PID。返回 nil 仅表示监控已成功建立；
// 实际等待由内核或后台 goroutine 完成。调用方负责在收到通知后执行优雅关闭或强制退出。
// 返回错误时调用方应终止启动，避免当前进程失去父进程生命周期约束。
func NotifyOnParentExit(parentPID int) (<-chan struct{}, error) {
	if parentPID <= 0 {
		err := fmt.Errorf("parent PID must be positive: %d", parentPID)
		console.Error("Failed to monitor parent process:", err)
		return nil, err
	}
	if err := checkParent(parentPID); err != nil {
		console.Error("Failed to monitor parent process:", err)
		return nil, err
	}
	parentExited, err := notifyOnParentExit(parentPID)
	if err != nil {
		console.Error("Failed to monitor parent process:", err)
		return nil, err
	}
	console.Debugf("Monitoring parent process %d", parentPID)
	return parentExited, nil
}

func checkParent(parentPID int) error {
	currentParentPID := os.Getppid()
	if currentParentPID != parentPID {
		return fmt.Errorf("%w: expected PID %d, current parent PID %d", ErrParentExited, parentPID, currentParentPID)
	}
	return nil
}
