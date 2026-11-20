//go:build linux

package process

import (
	"fmt"
	"os"
	"os/signal"
	"sync"
	"time"

	"golang.org/x/sys/unix"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// parentPollInterval 是进程级兜底监控的轮询间隔：pidfd 不可用时退化为
// 低频 PID/PPID 轮询，在退出响应速度和常驻唤醒开销之间取平衡。
const parentPollInterval = time.Second

func notifyOnParentExit(parentPID int) (<-chan struct{}, error) {
	parentExited := make(chan struct{})
	var closeOnce sync.Once
	notify := func() {
		closeOnce.Do(func() {
			console.Debugf("Parent process %d exited", parentPID)
			close(parentExited)
		})
	}

	// 快路径：PR_SET_PDEATHSIG 绑定的是调用 prctl 时的父进程线程。
	// 若 Python 从非主线程 spawn 且该线程先于进程退出，PDEATHSIG 只会触发一次
	// 伪通知，由下方的 PPID 复查过滤；此后不会再有信号，必须依赖进程级兜底。
	// SIGUSR1 保留给父死监控，signal.Notify 接管后不再触发默认终止行为，
	// 其他组件不得复用该信号。
	parentDeathSignals := make(chan os.Signal, 1)
	signal.Notify(parentDeathSignals, unix.SIGUSR1)
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
			// SIGUSR1 也可能来自其他进程或已退出的 spawn 线程；
			// 只有父 PID 改变才视为父进程退出。
			if checkParent(parentPID) == nil {
				continue
			}
			notify()
			return
		}
	}()

	// 进程级兜底：监控父进程本身而非创建本进程的线程。
	go func() {
		watchParentProcess(parentPID)
		notify()
	}()

	return parentExited, nil
}

// watchParentProcess 阻塞直到父进程退出。优先使用 pidfd（内核级通知），
// 不可用或中途出错时退化为低频 PPID 轮询。
func watchParentProcess(parentPID int) {
	if fd, err := unix.PidfdOpen(parentPID, 0); err == nil {
		defer func() { _ = unix.Close(fd) }()
		fds := []unix.PollFd{{Fd: int32(fd), Events: unix.POLLIN}}
		for {
			n, err := unix.Poll(fds, -1)
			if err == unix.EINTR {
				continue
			}
			if err != nil {
				// pidfd 等待异常（如被信号打断外的错误），退化为轮询兜底。
				break
			}
			if n > 0 {
				return // 进程退出，pidfd 可读
			}
		}
	}
	pollParentExit(parentPID)
}

// pollParentExit 以固定间隔检查 PPID 是否改变；父进程存活期间 PPID 保持不变，
// 整个进程退出后被重新托管（init/subreaper），PPID 随之改变。
// goroutine 启动前后存在极短竞态窗口，进入等待循环前先立即检查一次。
func pollParentExit(parentPID int) {
	if checkParent(parentPID) != nil {
		return
	}
	ticker := time.NewTicker(parentPollInterval)
	defer ticker.Stop()
	for range ticker.C {
		if checkParent(parentPID) != nil {
			return
		}
	}
}
