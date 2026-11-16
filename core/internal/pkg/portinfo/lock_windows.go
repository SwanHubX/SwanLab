//go:build windows

package portinfo

import (
	"errors"
	"fmt"
	"os"

	"golang.org/x/sys/windows"
)

const lockByteLength = 1

func openLockFile(lockPath string) (*os.File, error) {
	if linfo, err := os.Lstat(lockPath); err == nil {
		if linfo.Mode()&os.ModeSymlink != 0 || !linfo.Mode().IsRegular() {
			return nil, fmt.Errorf("lock path %s is not a regular file", lockPath)
		}
	} else if !errors.Is(err, os.ErrNotExist) {
		return nil, fmt.Errorf("lstat lock file %s: %w", lockPath, err)
	}

	f, err := os.OpenFile(lockPath, os.O_RDWR|os.O_CREATE, filePerm)
	if err != nil {
		return nil, fmt.Errorf("open lock file %s: %w", lockPath, err)
	}
	info, statErr := f.Stat()
	if statErr != nil {
		_ = f.Close()
		return nil, fmt.Errorf("stat lock file %s: %w", lockPath, statErr)
	}
	if !info.Mode().IsRegular() {
		_ = f.Close()
		return nil, fmt.Errorf("lock path %s is not a regular file", lockPath)
	}
	return f, nil
}

func tryFlockExclusive(f *os.File) (bool, error) {
	var ov windows.Overlapped
	err := windows.LockFileEx(
		windows.Handle(f.Fd()),
		windows.LOCKFILE_EXCLUSIVE_LOCK|windows.LOCKFILE_FAIL_IMMEDIATELY,
		0,
		lockByteLength,
		0,
		&ov,
	)
	if err == nil {
		return true, nil
	}
	if errors.Is(err, windows.ERROR_LOCK_VIOLATION) {
		return false, nil
	}
	return false, err
}

func unlockFile(f *os.File) error {
	var ov windows.Overlapped
	unlockErr := windows.UnlockFileEx(windows.Handle(f.Fd()), 0, lockByteLength, 0, &ov)
	closeErr := f.Close()
	if unlockErr != nil {
		return unlockErr
	}
	return closeErr
}
