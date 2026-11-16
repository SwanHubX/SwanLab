package portinfo

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"time"
)

// 启动锁默认值与文件格式常量。
const (
	LockFileName            = "core.lock"
	DefaultLockTimeout      = 35 * time.Second
	DefaultLockStaleTTL     = 30 * time.Second
	defaultLockPollInterval = 20 * time.Millisecond

	lockKeyPID       = "pid"
	lockKeyCreatedAt = "created_at"
	lockFieldCount   = 2
	decimalBase      = 10
)

// ErrLockTimeout 表示在超时时间内未能获取启动锁。
var ErrLockTimeout = errors.New("timed out waiting for core startup lock")

// LockOptions 配置启动锁的超时、陈旧判定与轮询行为；零值字段使用默认实现。
type LockOptions struct {
	Timeout        time.Duration
	StaleTTL       time.Duration
	PollInterval   time.Duration
	Now            func() time.Time
	IsProcessAlive func(int) bool
	Wait           func(context.Context, time.Duration) error
}

// FileLock 表示已持有的文件启动锁；Release 幂等且只删除归属本实例的锁文件。
type FileLock struct {
	path      string
	pid       int
	createdAt int64

	mu       sync.Mutex
	released bool
}

type lockInfo struct {
	pid       int
	createdAt int64
}

// Path 返回锁文件完整路径。
func (l *FileLock) Path() string {
	return l.path
}

// PID 返回锁记录的持有者进程 PID。
func (l *FileLock) PID() int {
	return l.pid
}

// Release 幂等释放启动锁；若锁已被他人接管则不误删。
func (l *FileLock) Release() error {
	l.mu.Lock()
	if l.released {
		l.mu.Unlock()
		return nil
	}
	l.released = true
	l.mu.Unlock()

	info, err := os.Lstat(l.path)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return nil
		}
		return fmt.Errorf("stat lock file %s: %w", l.path, err)
	}
	if info.Mode()&os.ModeSymlink != 0 || !info.Mode().IsRegular() {
		return fmt.Errorf("lock path %s is not a regular file", l.path)
	}

	current, err := parseLockFile(l.path)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return nil
		}
		return nil
	}
	if current.pid != l.pid || current.createdAt != l.createdAt {
		return nil
	}
	if err := os.Remove(l.path); err != nil && !errors.Is(err, os.ErrNotExist) {
		return fmt.Errorf("remove lock file %s: %w", l.path, err)
	}
	return nil
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

// TryAcquireLock 执行一次非阻塞抢锁；若已有存活持有者则返回 (nil, false, nil)。
func TryAcquireLock(runtimeDir string, opts LockOptions) (*FileLock, bool, error) {
	if err := ValidateRuntimeDir(runtimeDir); err != nil {
		return nil, false, err
	}
	cfg := normalizeLockOptions(opts)
	lockPath := filepath.Join(filepath.Clean(runtimeDir), LockFileName)
	pid := os.Getpid()
	now := cfg.Now()
	createdAt := now.UnixNano()

	err := createLockFile(lockPath, pid, createdAt)
	if err == nil {
		return &FileLock{path: lockPath, pid: pid, createdAt: createdAt}, true, nil
	}
	if !errors.Is(err, os.ErrExist) {
		return nil, false, fmt.Errorf("create lock file %s: %w", lockPath, err)
	}

	reclaimed, reclaimErr := reclaimIfStale(lockPath, now, cfg)
	if reclaimErr != nil {
		return nil, false, reclaimErr
	}
	if !reclaimed {
		return nil, false, nil
	}

	now = cfg.Now()
	createdAt = now.UnixNano()
	if retryErr := createLockFile(lockPath, pid, createdAt); retryErr != nil {
		if errors.Is(retryErr, os.ErrExist) {
			return nil, false, nil
		}
		return nil, false, fmt.Errorf("create lock file %s after reclaim: %w", lockPath, retryErr)
	}
	return &FileLock{path: lockPath, pid: pid, createdAt: createdAt}, true, nil
}

func normalizeLockOptions(opts LockOptions) LockOptions {
	if opts.Timeout <= 0 {
		opts.Timeout = DefaultLockTimeout
	}
	if opts.StaleTTL <= 0 {
		opts.StaleTTL = DefaultLockStaleTTL
	}
	if opts.PollInterval <= 0 {
		opts.PollInterval = defaultLockPollInterval
	}
	if opts.Now == nil {
		opts.Now = time.Now
	}
	if opts.IsProcessAlive == nil {
		opts.IsProcessAlive = processAlive
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

func createLockFile(lockPath string, pid int, createdAt int64) error {
	payload := marshalLock(pid, createdAt)
	dir := filepath.Dir(lockPath)
	tmp, err := os.CreateTemp(dir, ".core.lock-*")
	if err == nil {
		tmpName := tmp.Name()
		defer func() { _ = os.Remove(tmpName) }()

		writeErr := writeAll(tmp, payload)
		closeErr := tmp.Close()
		if writeErr != nil {
			return writeErr
		}
		if closeErr != nil {
			return closeErr
		}
		if chmodErr := os.Chmod(tmpName, filePerm); chmodErr != nil {
			return chmodErr
		}
		// 通过硬链接原子发布完整内容，避免先建空文件再写入期间被误判为损坏锁。
		linkErr := os.Link(tmpName, lockPath)
		if linkErr == nil || errors.Is(linkErr, os.ErrExist) {
			return linkErr
		}
	}

	// 少数不支持硬链接的文件系统回退到 O_CREATE|O_EXCL。
	f, err := os.OpenFile(lockPath, os.O_WRONLY|os.O_CREATE|os.O_EXCL, filePerm)
	if err != nil {
		return err
	}
	writeErr := writeAll(f, payload)
	closeErr := f.Close()
	if writeErr != nil {
		_ = os.Remove(lockPath)
		return writeErr
	}
	if closeErr != nil {
		_ = os.Remove(lockPath)
		return closeErr
	}
	return nil
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

func parseLockFile(path string) (lockInfo, error) {
	f, err := os.Open(path)
	if err != nil {
		return lockInfo{}, err
	}
	defer func() { _ = f.Close() }()

	data, err := io.ReadAll(io.LimitReader(f, MaxFileSize+1))
	if err != nil {
		return lockInfo{}, err
	}
	return parseLock(data)
}

func parseLock(data []byte) (lockInfo, error) {
	if len(data) > MaxFileSize {
		return lockInfo{}, fmt.Errorf("lock file exceeds %d bytes", MaxFileSize)
	}
	content := strings.TrimSuffix(string(data), "\n")
	if !strings.HasSuffix(content, "\n"+eofMarker) {
		return lockInfo{}, errors.New("lock file must end with an EOF line")
	}
	body := strings.TrimSuffix(content, "\n"+eofMarker)

	var parsed lockInfo
	seen := make(map[string]bool, lockFieldCount)
	for _, line := range strings.Split(body, "\n") {
		if line == "" {
			return lockInfo{}, errors.New("lock file contains empty line")
		}
		key, value, ok := strings.Cut(line, "=")
		if !ok {
			return lockInfo{}, fmt.Errorf("lock file line is not key=value: %q", line)
		}
		if seen[key] {
			return lockInfo{}, fmt.Errorf("lock file contains duplicate key %q", key)
		}
		switch key {
		case lockKeyPID:
			pid, err := parsePID(value)
			if err != nil {
				return lockInfo{}, err
			}
			parsed.pid = pid
		case lockKeyCreatedAt:
			ts, err := strconv.ParseInt(value, decimalBase, 64)
			if err != nil || ts <= 0 {
				return lockInfo{}, fmt.Errorf("invalid created_at %q", value)
			}
			parsed.createdAt = ts
		default:
			return lockInfo{}, fmt.Errorf("lock file contains unknown key %q", key)
		}
		seen[key] = true
	}
	if !seen[lockKeyPID] || !seen[lockKeyCreatedAt] {
		return lockInfo{}, errors.New("lock file missing required fields")
	}
	return parsed, nil
}

func reclaimIfStale(lockPath string, now time.Time, opts LockOptions) (bool, error) {
	info, err := os.Lstat(lockPath)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return true, nil
		}
		return false, fmt.Errorf("stat lock file %s: %w", lockPath, err)
	}
	if info.Mode()&os.ModeSymlink != 0 || !info.Mode().IsRegular() {
		return false, fmt.Errorf("lock path %s is not a regular file", lockPath)
	}

	parsed, parseErr := parseLockFile(lockPath)
	if parseErr != nil {
		if errors.Is(parseErr, os.ErrNotExist) {
			return true, nil
		}
		// 刚创建但尚未写完的锁文件不立即判废，只有超过 StaleTTL 才按崩溃残留清理。
		if now.Sub(info.ModTime()) < opts.StaleTTL {
			return false, nil
		}
	} else {
		created := time.Unix(0, parsed.createdAt)
		expired := now.Sub(created) >= opts.StaleTTL
		alive := opts.IsProcessAlive(parsed.pid)
		if !expired && alive {
			return false, nil
		}
	}

	latest, statErr := os.Lstat(lockPath)
	if statErr != nil {
		if errors.Is(statErr, os.ErrNotExist) {
			return true, nil
		}
		return false, fmt.Errorf("re-stat lock file %s: %w", lockPath, statErr)
	}
	if !os.SameFile(info, latest) {
		return false, nil
	}
	if rmErr := os.Remove(lockPath); rmErr != nil && !errors.Is(rmErr, os.ErrNotExist) {
		return false, fmt.Errorf("remove stale lock %s: %w", lockPath, rmErr)
	}
	return true, nil
}
