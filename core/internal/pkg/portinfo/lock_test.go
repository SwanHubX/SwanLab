package portinfo

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"runtime"
	"sync"
	"sync/atomic"
	"testing"
	"time"
)

func newTestRuntimeDir(t *testing.T) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), RuntimeDirPrefix+"test")
	if err := EnsureRuntimeDir(dir); err != nil {
		t.Fatalf("EnsureRuntimeDir: %v", err)
	}
	return dir
}

func TestTryAcquireAndReleaseLock(t *testing.T) {
	dir := newTestRuntimeDir(t)
	now := time.Date(2026, 11, 16, 14, 0, 0, 0, time.UTC)

	lock, acquired, err := TryAcquireLock(dir, LockOptions{Now: func() time.Time { return now }})
	if err != nil {
		t.Fatalf("TryAcquireLock: %v", err)
	}
	if !acquired || lock == nil {
		t.Fatal("expected lock to be acquired")
	}
	if lock.PID() != os.Getpid() {
		t.Fatalf("lock PID = %d, want %d", lock.PID(), os.Getpid())
	}
	if _, statErr := os.Stat(lock.Path()); statErr != nil {
		t.Fatalf("lock file missing after acquire: %v", statErr)
	}

	if relErr := lock.Release(); relErr != nil {
		t.Fatalf("first Release: %v", relErr)
	}
	if relErr := lock.Release(); relErr != nil {
		t.Fatalf("second Release: %v", relErr)
	}
	if _, statErr := os.Stat(lock.Path()); !errors.Is(statErr, os.ErrNotExist) {
		t.Fatalf("lock file still exists after release: %v", statErr)
	}
}

func TestReleaseDoesNotRemoveReplacedLock(t *testing.T) {
	dir := newTestRuntimeDir(t)
	base := time.Date(2026, 11, 16, 14, 0, 0, 0, time.UTC)

	oldLock, acquired, err := TryAcquireLock(dir, LockOptions{Now: func() time.Time { return base }})
	if err != nil || !acquired {
		t.Fatalf("first TryAcquireLock: acquired=%v err=%v", acquired, err)
	}

	// 模拟锁文件被新持有者接管覆盖。
	if rmErr := os.Remove(oldLock.Path()); rmErr != nil {
		t.Fatal(rmErr)
	}
	newTime := base.Add(time.Second)
	if createErr := createLockFile(oldLock.Path(), os.Getpid()+1, newTime.UnixNano()); createErr != nil {
		t.Fatal(createErr)
	}

	if relErr := oldLock.Release(); relErr != nil {
		t.Fatalf("Release on replaced lock: %v", relErr)
	}
	if _, statErr := os.Stat(oldLock.Path()); statErr != nil {
		t.Fatalf("new lock file was removed by old holder: %v", statErr)
	}
}

func TestTryAcquireRespectsLiveHolderAndReclaimsStale(t *testing.T) {
	dir := newTestRuntimeDir(t)
	now := time.Date(2026, 11, 16, 14, 0, 0, 0, time.UTC)
	clock := func() time.Time { return now }

	if err := createLockFile(filepath.Join(dir, LockFileName), 424242, now.UnixNano()); err != nil {
		t.Fatal(err)
	}

	// 持有者存活且未超时：不可抢占。
	lock, acquired, err := TryAcquireLock(dir, LockOptions{
		Now:            clock,
		StaleTTL:       30 * time.Second,
		IsProcessAlive: func(int) bool { return true },
	})
	if err != nil || acquired || lock != nil {
		t.Fatalf("expected busy lock, got lock=%v acquired=%v err=%v", lock, acquired, err)
	}

	// 持有者已死：立即回收并抢锁成功。
	lock, acquired, err = TryAcquireLock(dir, LockOptions{
		Now:            clock,
		StaleTTL:       30 * time.Second,
		IsProcessAlive: func(int) bool { return false },
	})
	if err != nil || !acquired || lock == nil {
		t.Fatalf("expected dead holder lock reclaim, got lock=%v acquired=%v err=%v", lock, acquired, err)
	}
	_ = lock.Release()

	// 持有者报告存活但超过 StaleTTL：超时回收。
	if createErr := createLockFile(filepath.Join(dir, LockFileName), 424242, now.UnixNano()); createErr != nil {
		t.Fatal(createErr)
	}
	now = now.Add(31 * time.Second)
	lock, acquired, err = TryAcquireLock(dir, LockOptions{
		Now:            clock,
		StaleTTL:       30 * time.Second,
		IsProcessAlive: func(int) bool { return true },
	})
	if err != nil || !acquired || lock == nil {
		t.Fatalf("expected expired lock reclaim, got lock=%v acquired=%v err=%v", lock, acquired, err)
	}
	_ = lock.Release()
}

func TestTryAcquireConcurrentSingleWinner(t *testing.T) {
	dir := newTestRuntimeDir(t)
	now := time.Now()

	const workers = 16
	start := make(chan struct{})
	var winners atomic.Int32
	var wg sync.WaitGroup
	wg.Add(workers)

	for range workers {
		go func() {
			defer wg.Done()
			<-start
			lock, acquired, err := TryAcquireLock(dir, LockOptions{
				Now:            func() time.Time { return now },
				IsProcessAlive: func(int) bool { return true },
			})
			if err != nil {
				t.Errorf("TryAcquireLock error: %v", err)
				return
			}
			if acquired {
				winners.Add(1)
				if lock == nil {
					t.Error("acquired=true returned nil lock")
				}
			}
		}()
	}

	close(start)
	wg.Wait()

	if got := winners.Load(); got != 1 {
		t.Fatalf("concurrent winners = %d, want 1", got)
	}
}

func TestAcquireLockWaitsAndSucceedsWhenReleased(t *testing.T) {
	dir := newTestRuntimeDir(t)
	now := time.Date(2026, 11, 16, 14, 0, 0, 0, time.UTC)
	lockPath := filepath.Join(dir, LockFileName)

	if err := createLockFile(lockPath, 999999, now.UnixNano()); err != nil {
		t.Fatal(err)
	}

	waits := 0
	lock, err := AcquireLock(context.Background(), dir, LockOptions{
		Timeout:        time.Second,
		PollInterval:   10 * time.Millisecond,
		Now:            func() time.Time { return now },
		IsProcessAlive: func(int) bool { return true },
		Wait: func(_ context.Context, d time.Duration) error {
			waits++
			now = now.Add(d)
			// 第一次等待时模拟前一个持有者释放锁。
			return os.Remove(lockPath)
		},
	})
	if err != nil {
		t.Fatalf("AcquireLock: %v", err)
	}
	defer func() { _ = lock.Release() }()
	if waits != 1 {
		t.Fatalf("waits = %d, want 1", waits)
	}
}

func TestAcquireLockTimesOut(t *testing.T) {
	dir := newTestRuntimeDir(t)
	now := time.Date(2026, 11, 16, 14, 0, 0, 0, time.UTC)

	if err := createLockFile(filepath.Join(dir, LockFileName), 999999, now.UnixNano()); err != nil {
		t.Fatal(err)
	}

	_, err := AcquireLock(context.Background(), dir, LockOptions{
		Timeout:        50 * time.Millisecond,
		StaleTTL:       time.Minute,
		PollInterval:   20 * time.Millisecond,
		Now:            func() time.Time { return now },
		IsProcessAlive: func(int) bool { return true },
		Wait: func(_ context.Context, d time.Duration) error {
			now = now.Add(d)
			return nil
		},
	})
	if !errors.Is(err, ErrLockTimeout) {
		t.Fatalf("err = %v, want ErrLockTimeout", err)
	}
}

func TestTryAcquireRejectsSymlinkLock(t *testing.T) {
	if runtime.GOOS == "windows" {
		t.Skip("symlink requires elevated privilege on Windows")
	}
	dir := newTestRuntimeDir(t)
	target := filepath.Join(t.TempDir(), "outside")
	if err := os.WriteFile(target, []byte("x"), 0o600); err != nil {
		t.Fatal(err)
	}
	if err := os.Symlink(target, filepath.Join(dir, LockFileName)); err != nil {
		t.Fatal(err)
	}

	_, acquired, err := TryAcquireLock(dir, LockOptions{})
	if err == nil || acquired {
		t.Fatalf("expected symlink lock rejection, got acquired=%v err=%v", acquired, err)
	}
}

func TestParseLockRejectsMalformed(t *testing.T) {
	cases := map[string]string{
		"missing EOF":     "pid=12\ncreated_at=100\n",
		"missing pid":     "created_at=100\nEOF\n",
		"missing created": "pid=12\nEOF\n",
		"bad pid":         "pid=0\ncreated_at=100\nEOF\n",
		"bad created_at":  "pid=12\ncreated_at=0\nEOF\n",
		"duplicate key":   "pid=12\npid=12\ncreated_at=100\nEOF\n",
		"unknown key":     "pid=12\ncreated_at=100\nextra=1\nEOF\n",
		"empty line":      "pid=12\n\ncreated_at=100\nEOF\n",
		"not key value":   "pid\ncreated_at=100\nEOF\n",
	}
	for name, content := range cases {
		if _, err := parseLock([]byte(content)); err == nil {
			t.Fatalf("%s: parseLock unexpectedly succeeded", name)
		}
	}
}
