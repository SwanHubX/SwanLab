// Package mime 按扩展名推断 MIME 类型，行为与 Python SDK 的 pkg.mime.guess_type 对齐。
//
// 自定义注册表覆盖 Python 显式 add_type 的扩展名，避免各平台 /etc/mime.types
// 差异导致同一份文件在两套 SDK 上得到不同的 Content-Type；其余扩展名交给标准库。
package mime

import (
	stdmime "mime"
	"path/filepath"
	"strings"
)

// DefaultType 是无法识别扩展名时的回退类型。
const DefaultType = "application/octet-stream"

// custom 与 Python pkg/mime/__init__.py 的 add_type 列表逐项对应。
var custom = map[string]string{
	".yaml":     "application/yaml",
	".yml":      "application/yaml",
	".png":      "image/png",
	".jpg":      "image/jpeg",
	".jpeg":     "image/jpeg",
	".gif":      "image/gif",
	".bmp":      "image/bmp",
	".webp":     "image/webp",
	".svg":      "image/svg+xml",
	".svgz":     "image/svg+xml",
	".tiff":     "image/tiff",
	".tif":      "image/tiff",
	".py":       "text/x-python",
	".txt":      "text/plain",
	".json":     "application/json",
	".md":       "text/markdown",
	".markdown": "text/markdown",
}

// Type 返回 path 的 MIME 类型；无法识别时回退 application/octet-stream。
// 扩展名大小写不敏感。
func Type(path string) string {
	ext := strings.ToLower(filepath.Ext(path))
	if mimeType, ok := custom[ext]; ok {
		return mimeType
	}
	if mimeType := stdmime.TypeByExtension(ext); mimeType != "" {
		return mimeType
	}
	return DefaultType
}
