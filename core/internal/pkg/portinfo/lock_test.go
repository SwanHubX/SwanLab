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
	now := time.Now()

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

	// 持锁期间二次抢锁必须失败。
	second, acquired2, err2 := TryAcquireLock(dir, LockOptions{})
	if err2 != nil || acquired2 || second != nil {
		t.Fatalf("expected busy lock while held, got lock=%v acquired=%v err=%v", second, acquired2, err2)
	}

	if relErr := lock.Release(); relErr != nil {
		t.Fatalf("first Release: %v", relErr)
	}
	if relErr := lock.Release(); relErr != nil {
		t.Fatalf("second Release: %v", relErr)
	}

	// 释放后再次抢锁成功，且重复调用旧锁 Release 不破坏新锁持有。
	next, acquired3, err3 := TryAcquireLock(dir, LockOptions{})
	if err3 != nil || !acquired3 || next == nil {
		t.Fatalf("re-acquire after release: lock=%v acquired=%v err=%v", next, acquired3, err3)
	}
	defer func() { _ = next.Release() }()
	if relErr := lock.Release(); relErr != nil {
		t.Fatalf("stale Release after takeover: %v", relErr)
	}

	third, acquired4, err4 := TryAcquireLock(dir, LockOptions{})
	if err4 != nil || acquired4 || third != nil {
		t.Fatalf("new lock was broken by old Release, got lock=%v acquired=%v err=%v", third, acquired4, err4)
	}
}

func TestTryAcquireConcurrentSingleWinner(t *testing.T) {
	dir := newTestRuntimeDir(t)

	const workers = 16
	start := make(chan struct{})
	var winners atomic.Int32
	var heldLock *FileLock
	var mu sync.Mutex
	var wg sync.WaitGroup
	wg.Add(workers)

	for range workers {
		go func() {
			defer wg.Done()
			<-start
			lock, acquired, err := TryAcquireLock(dir, LockOptions{})
			if err != nil {
				t.Errorf("TryAcquireLock error: %v", err)
				return
			}
			if acquired {
				winners.Add(1)
				mu.Lock()
				heldLock = lock
				mu.Unlock()
			}
		}()
	}

	close(start)
	wg.Wait()

	if heldLock != nil {
		defer func() { _ = heldLock.Release() }()
	}
	if got := winners.Load(); got != 1 {
		t.Fatalf("concurrent winners = %d, want 1", got)
	}
}

func TestLockInterleavedNoDoubleAcquire(t *testing.T) {
	dir := newTestRuntimeDir(t)
	first, acquired1, err1 := TryAcquireLock(dir, LockOptions{})
	if err1 != nil || !acquired1 {
		t.Fatalf("first TryAcquireLock: acquired=%v err=%v", acquired1, err1)
	}

	// 交错场景：A 仍在持锁与准备释放之间时，B 多次尝试均不得拿到锁；仅 A 真正 Release 后 B 才成功。
	bBlocked, bAcquired, bErr := TryAcquireLock(dir, LockOptions{})
	if bErr != nil || bAcquired || bBlocked != nil {
		t.Fatalf("interleaved B must fail while A holds lock, got acquired=%v err=%v", bAcquired, bErr)
	}

	if relErr := first.Release(); relErr != nil {
		t.Fatalf("first Release: %v", relErr)
	}

	second, acquired2, err2 := TryAcquireLock(dir, LockOptions{})
	if err2 != nil || !acquired2 || second == nil {
		t.Fatalf("second TryAcquireLock after release: acquired=%v err=%v", acquired2, err2)
	}
	defer func() { _ = second.Release() }()

	// A 再次幂等 Release 不影响 B 持锁。
	if relErr := first.Release(); relErr != nil {
		t.Fatalf("idempotent Release: %v", relErr)
	}
	cBlocked, cAcquired, cErr := TryAcquireLock(dir, LockOptions{})
	if cErr != nil || cAcquired || cBlocked != nil {
		t.Fatalf("third caller must not acquire while B holds lock, got acquired=%v err=%v", cAcquired, cErr)
	}
}

func TestAcquireLockWaitsAndSucceedsWhenReleased(t *testing.T) {
	dir := newTestRuntimeDir(t)
	holder, acquired, err := TryAcquireLock(dir, LockOptions{})
	if err != nil || !acquired {
		t.Fatalf("holder TryAcquireLock: acquired=%v err=%v", acquired, err)
	}

	now := time.Now()
	waits := 0
	lock, err := AcquireLock(context.Background(), dir, LockOptions{
		Timeout:      time.Second,
		PollInterval: 10 * time.Millisecond,
		Now:          func() time.Time { return now },
		Wait: func(_ context.Context, d time.Duration) error {
			waits++
			now = now.Add(d)
			return holder.Release()
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
	holder, acquired, err := TryAcquireLock(dir, LockOptions{})
	if err != nil || !acquired {
		t.Fatalf("holder TryAcquireLock: acquired=%v err=%v", acquired, err)
	}
	defer func() { _ = holder.Release() }()

	now := time.Now()
	_, err = AcquireLock(context.Background(), dir, LockOptions{
		Timeout:      50 * time.Millisecond,
		PollInterval: 20 * time.Millisecond,
		Now:          func() time.Time { return now },
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
