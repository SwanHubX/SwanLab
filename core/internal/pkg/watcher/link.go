package watcher

import (
	"context"
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

// copyChunk 是回退拷贝每次读取的块大小：块越小，取消响应越快。
const copyChunk = 1 << 20

// CreateSaveLinks 在 filesDir 下建立本地镜像并填充 TargetPath。
// 优先软链接，不支持时复制；已存在的目标不覆盖，返回新建数量。
// 复制回退前校验文件不超过 maxSize（<=0 表示不限制），拷贝过程响应 ctx 取消。
func CreateSaveLinks(ctx context.Context, saves []*savev1.SaveRecord, filesDir string, maxSize int64) (int, error) {
	if filesDir == "" {
		return 0, errors.New("watcher: files dir is required")
	}
	abs, err := filepath.Abs(filesDir)
	if err != nil {
		return 0, fmt.Errorf("watcher: resolve files dir: %w", err)
	}
	created := 0
	for _, save := range saves {
		if err = ctx.Err(); err != nil {
			return created, err
		}
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
		if err = linkFile(ctx, save.GetSourcePath(), target, maxSize); err != nil {
			return created, err
		}
		created++
	}
	return created, nil
}

// linkFile 优先建软链接；文件系统不支持时回退为整文件拷贝。
func linkFile(ctx context.Context, source, target string, maxSize int64) error {
	if err := os.Symlink(source, target); err == nil {
		return nil
	} else if copyErr := copyFile(ctx, source, target, maxSize); copyErr != nil {
		return fmt.Errorf("link %q -> %q: symlink: %v; copy: %w", target, source, err, copyErr)
	}
	return nil
}

// copyFile 整文件拷贝，按块响应取消；超过 maxSize 的文件直接拒绝，避免白拷贝一遍再被上传端丢弃。
func copyFile(ctx context.Context, source, target string, maxSize int64) error {
	in, err := os.Open(source)
	if err != nil {
		return err
	}
	defer func() { _ = in.Close() }()

	info, err := in.Stat()
	if err != nil {
		return err
	}
	if maxSize > 0 && info.Size() > maxSize {
		return fmt.Errorf("%q is %d bytes, exceeds the save size limit of %d bytes", source, info.Size(), maxSize)
	}
	out, err := os.OpenFile(target, os.O_WRONLY|os.O_CREATE|os.O_EXCL, info.Mode().Perm())
	if err != nil {
		return err
	}
	if _, err = io.CopyBuffer(out, ctxReader{ctx: ctx, reader: in}, make([]byte, copyChunk)); err != nil {
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

// ctxReader 在每次读取前检查取消，让大文件拷贝能及时退出。
type ctxReader struct {
	ctx    context.Context
	reader io.Reader
}

// Read 实现 io.Reader。
func (r ctxReader) Read(p []byte) (int, error) {
	if err := r.ctx.Err(); err != nil {
		return 0, err
	}
	return r.reader.Read(p)
}
