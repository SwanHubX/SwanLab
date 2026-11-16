package portinfo

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"time"
)

// 启动锁默认值与文件格式常量。
const (
	LockFileName            = "core.lock"
	DefaultLockTimeout      = 35 * time.Second
	defaultLockPollInterval = 20 * time.Millisecond

	lockKeyPID       = "pid"
	lockKeyCreatedAt = "created_at"
)

// ErrLockTimeout 表示在超时时间内未能获取启动锁。
var ErrLockTimeout = errors.New("timed out waiting for core startup lock")

// LockOptions 配置启动锁的超时与轮询行为；零值字段使用默认实现。
type LockOptions struct {
	Timeout      time.Duration
	PollInterval time.Duration
	Now          func() time.Time
	Wait         func(context.Context, time.Duration) error
}

// FileLock 表示通过内核独占文件锁持有的启动锁；Release 幂等释放但不删除锁文件，避免 ABA 竞态。
type FileLock struct {
	path      string
	pid       int
	createdAt int64

	mu       sync.Mutex
	file     *os.File
	released bool
}

// Path 返回锁文件完整路径。
func (l *FileLock) Path() string {
	return l.path
}

// PID 返回锁记录的持有者进程 PID。
func (l *FileLock) PID() int {
	return l.pid
}

// Release 幂等解锁并关闭底层文件句柄；不删除锁文件本身。
func (l *FileLock) Release() error {
	l.mu.Lock()
	if l.released {
		l.mu.Unlock()
		return nil
	}
	l.released = true
	f := l.file
	l.file = nil
	l.mu.Unlock()

	if f == nil {
		return nil
	}
	return unlockFile(f)
}

// AcquireLock 在运行时目录内竞争启动锁，直到获取成功、上下文取消或超时。
func AcquireLock(ctx context.Context, runtimeDir string, opts LockOptions) (*FileLock, error) {
	cfg := normalizeLockOptions(opts)
	deadline := cfg.Now().Add(cfg.Timeout)

	for {
		if ctxErr := ctx.Err(); ctxErr != nil {
			return nil, fmt.Errorf("%w: %w", ErrLockTimeout, ctxErr)
		}

		lock, acquired, err := TryAcquireLock(runtimeDir, cfg)
		if err != nil {
			return nil, err
		}
		if acquired {
			return lock, nil
		}

		remaining := deadline.Sub(cfg.Now())
		if remaining <= 0 {
			return nil, ErrLockTimeout
		}
		waitDur := cfg.PollInterval
		if remaining < waitDur {
			waitDur = remaining
		}
		if waitErr := cfg.Wait(ctx, waitDur); waitErr != nil {
			return nil, fmt.Errorf("%w: %w", ErrLockTimeout, waitErr)
		}
	}
}

// TryAcquireLock 执行一次非阻塞内核独占加锁；若已被其它持有者占用返回 (nil, false, nil)。
func TryAcquireLock(runtimeDir string, opts LockOptions) (*FileLock, bool, error) {
	if err := ValidateRuntimeDir(runtimeDir); err != nil {
		return nil, false, err
	}
	cfg := normalizeLockOptions(opts)
	lockPath := filepath.Join(filepath.Clean(runtimeDir), LockFileName)

	f, err := openLockFile(lockPath)
	if err != nil {
		return nil, false, err
	}

	locked, lockErr := tryFlockExclusive(f)
	if lockErr != nil {
		_ = f.Close()
		return nil, false, fmt.Errorf("lock %s: %w", lockPath, lockErr)
	}
	if !locked {
		_ = f.Close()
		return nil, false, nil
	}

	pid := os.Getpid()
	createdAt := cfg.Now().UnixNano()
	if writeErr := writeLockMetadata(f, pid, createdAt); writeErr != nil {
		_ = unlockFile(f)
		return nil, false, fmt.Errorf("write lock metadata %s: %w", lockPath, writeErr)
	}

	return &FileLock{
		path:      lockPath,
		pid:       pid,
		createdAt: createdAt,
		file:      f,
	}, true, nil
}

func normalizeLockOptions(opts LockOptions) LockOptions {
	if opts.Timeout <= 0 {
		opts.Timeout = DefaultLockTimeout
	}
	if opts.PollInterval <= 0 {
		opts.PollInterval = defaultLockPollInterval
	}
	if opts.Now == nil {
		opts.Now = time.Now
	}
	if opts.Wait == nil {
		opts.Wait = waitWithContext
	}
	return opts
}

func waitWithContext(ctx context.Context, d time.Duration) error {
	timer := time.NewTimer(d)
	defer timer.Stop()
	select {
	case <-ctx.Done():
		return ctx.Err()
	case <-timer.C:
		return nil
	}
}

func writeLockMetadata(f *os.File, pid int, createdAt int64) error {
	if err := f.Truncate(0); err != nil {
		return err
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		return err
	}
	return writeAll(f, marshalLock(pid, createdAt))
}

func writeAll(f *os.File, data []byte) error {
	for len(data) > 0 {
		n, err := f.Write(data)
		if err != nil {
			return err
		}
		data = data[n:]
	}
	return nil
}

func marshalLock(pid int, createdAt int64) []byte {
	var b strings.Builder
	fmt.Fprintf(&b, "%s=%d\n", lockKeyPID, pid)
	fmt.Fprintf(&b, "%s=%d\n", lockKeyCreatedAt, createdAt)
	b.WriteString(eofMarker + "\n")
	return []byte(b.String())
}
