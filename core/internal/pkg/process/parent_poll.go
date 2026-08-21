//go:build !linux && !windows

package process

import (
	"os"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// parentCheckInterval 在退出响应速度和常驻唤醒开销之间取平衡。
const parentCheckInterval = time.Second

func notifyOnParentExit(parentPID int) (<-chan struct{}, error) {
	parentExited := make(chan struct{})

	// 这些平台没有统一可用的父进程死亡通知，因此使用低频后台轮询作为兜底。
	go func() {
		ticker := time.NewTicker(parentCheckInterval)
		defer ticker.Stop()

		waitForParentExit(parentPID, ticker.C, os.Getppid)
		console.Debugf("Parent process %d exited", parentPID)
		close(parentExited)
	}()
	return parentExited, nil
}

func waitForParentExit(parentPID int, ticks <-chan time.Time, getParentPID func() int) {
	// 启动 goroutine 前后仍存在极短竞态窗口，进入等待循环前先立即检查一次。
	if getParentPID() != parentPID {
		return
	}
	for range ticks {
		if getParentPID() != parentPID {
			return
		}
	}
}
