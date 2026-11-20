package watcher

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/swanhubx/swanlab/core/internal/pkg/fs"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

// dirPerm 是本地镜像目录的权限。
const dirPerm = 0o755

// CreateSaveLinks 在 filesDir 下建立本地镜像并填充 TargetPath。
// 优先软链接，不支持时复制；已存在的目标不覆盖，返回新建数量。
func CreateSaveLinks(saves []*savev1.SaveRecord, filesDir string) (int, error) {
	if filesDir == "" {
		return 0, errors.New("watcher: files dir is required")
	}
	abs, err := filepath.Abs(filesDir)
	if err != nil {
		return 0, fmt.Errorf("watcher: resolve files dir: %w", err)
	}
	created := 0
	for _, save := range saves {
		if save == nil {
			return created, errors.New("watcher: nil save record")
		}
		target := fs.MirrorPath(abs, save.GetName())
		save.TargetPath = target
		if _, err = os.Lstat(target); err == nil {
			continue
		} else if !errors.Is(err, os.ErrNotExist) {
			return created, fmt.Errorf("inspect %q: %w", target, err)
		}
		if err = os.MkdirAll(filepath.Dir(target), dirPerm); err != nil {
			return created, fmt.Errorf("create dir for %q: %w", target, err)
		}
		if err = linkFile(save.GetSourcePath(), target); err != nil {
			return created, err
		}
		created++
	}
	return created, nil
}

func linkFile(source, target string) error {
	if err := os.Symlink(source, target); err == nil {
		return nil
	} else if copyErr := copyFile(source, target); copyErr != nil {
		return fmt.Errorf("link %q -> %q: symlink: %v; copy: %w", target, source, err, copyErr)
	}
	return nil
}

func copyFile(source, target string) error {
	in, err := os.Open(source)
	if err != nil {
		return err
	}
	defer func() { _ = in.Close() }()

	info, err := in.Stat()
	if err != nil {
		return err
	}
	out, err := os.OpenFile(target, os.O_WRONLY|os.O_CREATE|os.O_EXCL, info.Mode().Perm())
	if err != nil {
		return err
	}
	if _, err = io.Copy(out, in); err != nil {
		_ = out.Close()
		_ = os.Remove(target)
		return err
	}
	if err = out.Close(); err != nil {
		_ = os.Remove(target)
		return err
	}
	return nil
}
