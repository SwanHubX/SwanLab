package portinfo

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// per-user 运行时目录与内部约定文件名。
const (
	RuntimeDirPrefix   = "swanlab-core-"
	CorePortFileName   = "core.port"
	CoreSocketFileName = "core.sock"
	CoreStderrFileName = "core.stderr.log"

	// RuntimeDirPerm 限定运行时目录仅当前用户可读写进入。
	RuntimeDirPerm os.FileMode = 0o700
)

// MaxUnixPathLen 是 POSIX UDS 路径字节上限，与 portinfo 解析校验保持一致。
const MaxUnixPathLen = maxUnixPathLen

// DefaultRuntimeDir 返回当前用户的固定运行时目录路径（不自动创建）。
func DefaultRuntimeDir() (string, error) {
	return runtimeDirInRoot(defaultRuntimeRoot())
}

func runtimeDirInRoot(root string) (string, error) {
	if strings.TrimSpace(root) == "" {
		return "", errors.New("runtime root directory is empty")
	}
	userSegment, err := currentUserSegment()
	if err != nil {
		return "", err
	}
	return filepath.Join(filepath.Clean(root), RuntimeDirPrefix+userSegment), nil
}

// EnsureRuntimeDir 确保运行时目录存在且满足当前用户的安全约束。
func EnsureRuntimeDir(dir string) error {
	if strings.TrimSpace(dir) == "" {
		return errors.New("runtime directory is empty")
	}
	cleaned := filepath.Clean(dir)
	err := os.Mkdir(cleaned, RuntimeDirPerm)
	if err == nil {
		// 显式 chmod 消除进程 umask 影响，确保新建目录严格为 0700。
		if chmodErr := ensureRuntimeDirPerm(cleaned); chmodErr != nil {
			return chmodErr
		}
	} else if !errors.Is(err, os.ErrExist) {
		return fmt.Errorf("create runtime dir %s: %w", cleaned, err)
	}
	return ValidateRuntimeDir(cleaned)
}

// ValidateRuntimeDir 校验运行时目录是当前用户拥有的真实目录（非符号链接），且 POSIX 下权限为 0700。
func ValidateRuntimeDir(dir string) error {
	if strings.TrimSpace(dir) == "" {
		return errors.New("runtime directory is empty")
	}
	cleaned := filepath.Clean(dir)
	info, err := os.Lstat(cleaned)
	if err != nil {
		return fmt.Errorf("stat runtime dir %s: %w", cleaned, err)
	}
	if info.Mode()&os.ModeSymlink != 0 {
		return fmt.Errorf("runtime dir %s must not be a symlink", cleaned)
	}
	if !info.IsDir() {
		return fmt.Errorf("runtime path %s is not a directory", cleaned)
	}
	return validateRuntimeDirOwner(cleaned, info)
}

// ValidateSocketPath 校验 POSIX 下默认 UDS 路径不超过 sun_path 长度限制。
func ValidateSocketPath(dir string) error {
	return validateSocketPathPlatform(filepath.Clean(dir))
}
