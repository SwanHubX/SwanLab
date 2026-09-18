// Package portinfo 实现 swanlab-core 端口信息文件（port-file）的严格读写。
//
// port-file 是 Go core 服务在 listen 成功后向 Python SDK 回报监听端点的约定文件，
// 同时携带 RPC 鉴权所需的 auth token。文件格式（v1）为若干行 key=value 文本，
// 以独立的 EOF 行结尾：
//
//	protocol=1
//	unix=/short/private/runtime/core.sock   （POSIX 平台）
//	sock=12345                              （Windows 平台，替代 unix 行）
//	auth=<base64url 编码的 256-bit 随机 token>
//	EOF
//
// 写入通过同目录临时文件 + fsync + chmod(0600) + rename 原子提交，读者只会看到
// 完整文件；解析对重复 key、未知协议版本、缺失字段、非法端口、非 EOF 结尾和
// 超长内容一律拒绝。auth token 属于敏感信息，任何错误消息中不得包含其值。
package portinfo

import (
	"crypto/rand"
	"encoding/base64"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// 协议与格式约束。
const (
	// ProtocolVersion 是当前 port-file 格式版本。
	ProtocolVersion = 1
	// MaxFileSize 是 port-file 内容长度上限，超出视为损坏或恶意构造。
	MaxFileSize = 4096
	// authTokenBytes 是 auth token 解码后的字节数（256-bit）。
	authTokenBytes = 32
	// maxUnixPathLen 限定 unix 端点路径长度，实际可用长度还受 sun_path 限制。
	maxUnixPathLen = 256
	// filePerm 是 port-file 的 owner-only 权限。
	filePerm = 0o600
	// eofMarker 是文件结尾标记，必须独占一行。
	eofMarker = "EOF"
)

// 字段名约定，解析与序列化共用。
const (
	keyProtocol = "protocol"
	keyUnix     = "unix"
	keySock     = "sock"
	keyAuth     = "auth"
)

// Info 是 port-file 的结构化内容，UnixPath 与 SockPort 二选一。
type Info struct {
	Protocol  int
	UnixPath  string
	SockPort  int
	AuthToken string
}

// NewAuthToken 生成 base64url 编码的 256-bit 随机 token，
// 用于写入 port-file 的 auth 字段。
func NewAuthToken() (string, error) {
	buf := make([]byte, authTokenBytes)
	if _, err := rand.Read(buf); err != nil {
		return "", fmt.Errorf("generate auth token: %w", err)
	}
	return base64.RawURLEncoding.EncodeToString(buf), nil
}

// Marshal 将 Info 序列化为 v1 格式字节串，序列化前完成全部校验。
func Marshal(info *Info) ([]byte, error) {
	if err := validate(info); err != nil {
		return nil, err
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s=%d\n", keyProtocol, info.Protocol)
	if info.UnixPath != "" {
		fmt.Fprintf(&b, "%s=%s\n", keyUnix, info.UnixPath)
	} else {
		fmt.Fprintf(&b, "%s=%d\n", keySock, info.SockPort)
	}
	fmt.Fprintf(&b, "%s=%s\n", keyAuth, info.AuthToken)
	b.WriteString(eofMarker + "\n")
	return []byte(b.String()), nil
}

// WriteFile 原子写入 port-file：先写同目录临时文件，fsync、chmod 后 rename 覆盖目标。
// 任一步失败都不会破坏既有目标文件，临时文件会被清理。
func WriteFile(path string, info *Info) error {
	data, err := Marshal(info)
	if err != nil {
		return err
	}
	dir := filepath.Dir(path)
	tmp, err := os.CreateTemp(dir, ".portinfo-*")
	if err != nil {
		return fmt.Errorf("create temp port-file in %s: %w", dir, err)
	}
	tmpName := tmp.Name()
	renamed := false
	defer func() {
		if !renamed {
			_ = os.Remove(tmpName)
		}
	}()
	if _, err := tmp.Write(data); err != nil {
		_ = tmp.Close()
		return fmt.Errorf("write temp port-file: %w", err)
	}
	// fsync 尽力而为：部分网络文件系统（如 smbfs）不支持 fsync 而返回 EINVAL。
	// port-file 为短命文件且由同机读者在秒级内读取，可见性一致性由 rename 原子性保证，
	// fsync 仅作为崩溃防护，失败不阻塞端点回报。
	_ = tmp.Sync()
	if err := tmp.Close(); err != nil {
		return fmt.Errorf("close temp port-file: %w", err)
	}
	if err := os.Chmod(tmpName, filePerm); err != nil {
		return fmt.Errorf("chmod temp port-file: %w", err)
	}
	if err := os.Rename(tmpName, path); err != nil {
		return fmt.Errorf("commit port-file: %w", err)
	}
	renamed = true
	syncDir(dir)
	return nil
}

// ParseFile 读取并严格解析 port-file。文件超过 MaxFileSize 时按损坏处理。
func ParseFile(path string) (Info, error) {
	f, err := os.Open(path)
	if err != nil {
		return Info{}, fmt.Errorf("open port-file: %w", err)
	}
	defer func() { _ = f.Close() }()
	data, err := io.ReadAll(io.LimitReader(f, MaxFileSize+1))
	if err != nil {
		return Info{}, fmt.Errorf("read port-file: %w", err)
	}
	return Parse(data)
}

// Parse 严格解析 v1 格式内容。键值行顺序不限，但每个必要字段恰好出现一次，
// 必须以独立 EOF 行结尾（末尾换行可选）。
func Parse(data []byte) (Info, error) {
	if len(data) > MaxFileSize {
		return Info{}, fmt.Errorf("port-file exceeds %d bytes", MaxFileSize)
	}
	content := strings.TrimSuffix(string(data), "\n")
	if !strings.HasSuffix(content, "\n"+eofMarker) && content != eofMarker {
		return Info{}, errors.New("port-file must end with an EOF line")
	}
	body := strings.TrimSuffix(content, "\n"+eofMarker)
	if content == eofMarker {
		body = ""
	}

	var info Info
	seen := make(map[string]bool, 3)
	for _, line := range strings.Split(body, "\n") {
		if line == "" {
			return Info{}, errors.New("port-file contains empty line")
		}
		key, value, ok := strings.Cut(line, "=")
		if !ok {
			return Info{}, fmt.Errorf("port-file line is not key=value: %q", line)
		}
		if seen[key] {
			return Info{}, fmt.Errorf("port-file contains duplicate key %q", key)
		}
		switch key {
		case keyProtocol:
			if value != strconv.Itoa(ProtocolVersion) {
				return Info{}, fmt.Errorf("unsupported port-file protocol %q", value)
			}
			info.Protocol = ProtocolVersion
		case keyUnix:
			if !strings.HasPrefix(value, "/") || len(value) > maxUnixPathLen {
				return Info{}, errors.New("unix endpoint must be an absolute path within length limit")
			}
			info.UnixPath = value
		case keySock:
			port, err := parsePort(value)
			if err != nil {
				return Info{}, err
			}
			info.SockPort = port
		case keyAuth:
			if err := validateAuthToken(value); err != nil {
				return Info{}, err
			}
			info.AuthToken = value
		default:
			return Info{}, fmt.Errorf("port-file contains unknown key %q", key)
		}
		seen[key] = true
	}

	if !seen[keyProtocol] {
		return Info{}, fmt.Errorf("port-file missing %q field", keyProtocol)
	}
	if !seen[keyAuth] {
		return Info{}, fmt.Errorf("port-file missing %q field", keyAuth)
	}
	if seen[keyUnix] == seen[keySock] {
		return Info{}, fmt.Errorf("port-file must contain exactly one of %q or %q", keyUnix, keySock)
	}
	return info, nil
}

// validate 校验 Info 的全部字段约束。
func validate(info *Info) error {
	if info.Protocol != ProtocolVersion {
		return fmt.Errorf("unsupported port-file protocol %d", info.Protocol)
	}
	hasUnix := info.UnixPath != ""
	hasSock := info.SockPort != 0
	if hasUnix == hasSock {
		return fmt.Errorf("exactly one of unix path or sock port must be set")
	}
	if hasUnix {
		if !strings.HasPrefix(info.UnixPath, "/") || len(info.UnixPath) > maxUnixPathLen {
			return fmt.Errorf("unix socket path must be absolute and at most %d bytes", maxUnixPathLen)
		}
	} else if info.SockPort < 1 || info.SockPort > 65535 {
		return fmt.Errorf("sock port %d out of range", info.SockPort)
	}
	return validateAuthToken(info.AuthToken)
}

// validateAuthToken 校验 token 是 base64url 编码且解码后恰好 256-bit。
// 错误消息不回显 token 值，避免敏感信息泄露。
func validateAuthToken(token string) error {
	raw, err := base64.RawURLEncoding.Strict().DecodeString(token)
	if err != nil || len(raw) != authTokenBytes {
		return errors.New("auth token must be base64url-encoded 256-bit value")
	}
	return nil
}

// parsePort 严格解析十进制端口号：仅数字、无前导零、范围 1-65535。
func parsePort(value string) (int, error) {
	if value == "" || strings.HasPrefix(value, "0") || !isDigits(value) {
		return 0, fmt.Errorf("invalid sock port %q", value)
	}
	port, err := strconv.Atoi(value)
	if err != nil || port < 1 || port > 65535 {
		return 0, fmt.Errorf("invalid sock port %q", value)
	}
	return port, nil
}

func isDigits(s string) bool {
	for _, r := range s {
		if r < '0' || r > '9' {
			return false
		}
	}
	return true
}

// syncDir 尽力持久化目录项，使 rename 结果落盘；失败不影响写入结果。
func syncDir(dir string) {
	d, err := os.Open(dir)
	if err != nil {
		return
	}
	defer func() { _ = d.Close() }()
	_ = d.Sync()
}
