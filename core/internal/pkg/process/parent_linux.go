//go:build linux

package process

import (
	"fmt"
	"os"
	"os/signal"

	"golang.org/x/sys/unix"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

func notifyOnParentExit(parentPID int) (<-chan struct{}, error) {
	parentExited := make(chan struct{})
	parentDeathSignals := make(chan os.Signal, 1)
	signal.Notify(parentDeathSignals, unix.SIGUSR1)

	// PR_SET_PDEATHSIG 绑定的是调用 prctl 时的父进程，而不是任意指定 PID。
	// 外层和下方的 PID 校验共同保证该父进程就是调用方预期的 parentPID。
	if err := unix.Prctl(unix.PR_SET_PDEATHSIG, uintptr(unix.SIGUSR1), 0, 0, 0); err != nil {
		signal.Stop(parentDeathSignals)
		return nil, fmt.Errorf("set parent death signal: %w", err)
	}

	// 父进程可能恰好在外层首次检查与 prctl 调用之间退出，必须再次确认。
	if err := checkParent(parentPID); err != nil {
		signal.Stop(parentDeathSignals)
		return nil, err
	}

	go func() {
		defer signal.Stop(parentDeathSignals)
		for range parentDeathSignals {
			// SIGUSR1 也可能来自其他进程；只有父 PID 改变才视为父进程退出。
			if checkParent(parentPID) == nil {
				continue
			}
			console.Debugf("Parent process %d exited", parentPID)
			close(parentExited)
			return
		}
	}()

	return parentExited, nil
}
