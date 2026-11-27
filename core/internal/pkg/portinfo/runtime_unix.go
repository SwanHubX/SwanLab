//go:build !windows

package portinfo

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"syscall"
)

const posixRuntimeRoot = "/tmp"

func defaultRuntimeRoot() string {
	return posixRuntimeRoot
}

func currentUserSegment() (string, error) {
	return strconv.Itoa(os.Getuid()), nil
}

func ensureRuntimeDirPerm(dir string) error {
	if err := os.Chmod(dir, RuntimeDirPerm); err != nil {
		return fmt.Errorf("chmod runtime dir %s: %w", dir, err)
	}
	return nil
}

func validateRuntimeDirOwner(dir string, info os.FileInfo) error {
	if perm := info.Mode().Perm(); perm != RuntimeDirPerm {
		return fmt.Errorf("runtime dir %s mode is %#o, want %#o", dir, perm, RuntimeDirPerm)
	}
	stat, ok := info.Sys().(*syscall.Stat_t)
	if !ok {
		return fmt.Errorf("cannot inspect owner of runtime dir %s", dir)
	}
	uid := os.Getuid()
	if int(stat.Uid) != uid {
		return fmt.Errorf("runtime dir %s is owned by uid %d, want %d", dir, stat.Uid, uid)
	}
	return nil
}

func validateSocketPathPlatform(dir string) error {
	sockPath := filepath.Join(dir, CoreSocketFileName)
	if len(sockPath) > maxUnixPathLen {
		return fmt.Errorf("unix socket path %s exceeds %d bytes", sockPath, maxUnixPathLen)
	}
	return nil
}
