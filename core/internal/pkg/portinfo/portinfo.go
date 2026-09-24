// Package portinfo 实现 swanlab-core 端口信息文件（port-file）的读写。
//
// port-file 是 core 服务在 listen 成功后向调用方回报监听端点的约定文件，
// 承担临时 endpoint discovery。文件格式（v1）为若干行 key=value 文本，
// 以独立的 EOF 行结尾：
//
//	protocol=1
//	pid=<core-pid>                        （诊断与 stale 检查辅助）
//	unix=/short/private/runtime/core.sock （POSIX 平台）
//	sock=12345                            （Windows 或 UDS 回退，替代 unix 行）
//	EOF
//
// 写入经同目录临时文件 + fsync + chmod(0600) + rename 原子提交；解析拒绝
// 重复 key、未知 key、未知协议版本、缺失字段、非法端口或 pid、非 EOF 结尾
// 和超长内容。
package portinfo

import (
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
	// maxPID 是 pid 字段的取值上限。
	maxPID = 1<<31 - 1
	// maxUnixPathLen 限定 unix 端点路径长度，实际可用长度还受 sun_path 限制。
	maxUnixPathLen = 256
	// filePerm 是 port-file 的 owner-only 权限。
	filePerm = 0o600
	// eofMarker 是文件结尾标记，独占一行。
	eofMarker = "EOF"
)

// 字段名约定，解析与序列化共用。
const (
	keyProtocol = "protocol"
	keyUnix     = "unix"
	keySock     = "sock"
	keyPID      = "pid"
)

// Info 是 port-file 的结构化内容，UnixPath 与 SockPort 二选一。
type Info struct {
	Protocol int
	UnixPath string
	SockPort int
	PID      int
}

// Marshal 将 Info 序列化为 v1 格式字节串，序列化前完成校验。
func Marshal(info *Info) ([]byte, error) {
	if err := validate(info); err != nil {
		return nil, err
	}
	var b strings.Builder
	fmt.Fprintf(&b, "%s=%d\n", keyProtocol, info.Protocol)
	fmt.Fprintf(&b, "%s=%d\n", keyPID, info.PID)
	if info.UnixPath != "" {
		fmt.Fprintf(&b, "%s=%s\n", keyUnix, info.UnixPath)
	} else {
		fmt.Fprintf(&b, "%s=%d\n", keySock, info.SockPort)
	}
	b.WriteString(eofMarker + "\n")
	return []byte(b.String()), nil
}

// WriteFile 原子写入 port-file：先写同目录临时文件，fsync、chmod 后 rename 覆盖目标。
// 任一步失败不破坏既有目标文件，临时文件会被清理。
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

// ParseFile 读取并解析 port-file。文件超过 MaxFileSize 时按损坏处理。
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

// Parse 解析 v1 格式内容。键值行顺序不限，每个字段出现一次，
// 以独立 EOF 行结尾（末尾换行可选）。
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
	seen := make(map[string]bool, 4)
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
		case keyPID:
			pid, err := parsePID(value)
			if err != nil {
				return Info{}, err
			}
			info.PID = pid
		default:
			return Info{}, fmt.Errorf("port-file contains unknown key %q", key)
		}
		seen[key] = true
	}

	if !seen[keyProtocol] {
		return Info{}, fmt.Errorf("port-file missing %q field", keyProtocol)
	}
	if !seen[keyPID] {
		return Info{}, fmt.Errorf("port-file missing %q field", keyPID)
	}
	if seen[keyUnix] == seen[keySock] {
		return Info{}, fmt.Errorf("port-file must contain exactly one of %q or %q", keyUnix, keySock)
	}
	return info, nil
}

// validate 校验 Info 的字段约束。
func validate(info *Info) error {
	if info.Protocol != ProtocolVersion {
		return fmt.Errorf("unsupported port-file protocol %d", info.Protocol)
	}
	if info.PID < 1 || info.PID > maxPID {
		return fmt.Errorf("pid %d out of range", info.PID)
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
	return nil
}

// parsePID 解析十进制 pid：数字、无前导零、范围 1-maxPID。
func parsePID(value string) (int, error) {
	if value == "" || strings.HasPrefix(value, "0") || !isDigits(value) {
		return 0, fmt.Errorf("invalid pid %q", value)
	}
	pid, err := strconv.Atoi(value)
	if err != nil || pid < 1 || pid > maxPID {
		return 0, fmt.Errorf("invalid pid %q", value)
	}
	return pid, nil
}

// parsePort 解析十进制端口号：数字、无前导零、范围 1-65535。
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

// syncDir 持久化目录项，使 rename 结果落盘；失败不影响写入结果。
func syncDir(dir string) {
	d, err := os.Open(dir)
	if err != nil {
		return
	}
	defer func() { _ = d.Close() }()
	_ = d.Sync()
}
