// Package fs 提供文件名安全化工具，把任意名称映射为文件系统兼容的形式。
// 行为与 Python SDK 的 fs.safe_fmt / fs.safe_truncate 保持一致。
package fs

import (
	"fmt"
	"strings"
)

const (
	truncHead     = 3
	truncEllipsis = "..."
	// minTruncateLen 是 SafeTruncate 允许的最小 maxLen。
	minTruncateLen = truncHead + len(truncEllipsis) + 1
)

// osUnsafe 判断文件名中不允许的字符：<>:"/\|?*#% 以及 C0 控制字符与 DEL。
func osUnsafe(r rune) bool {
	switch r {
	case '<', '>', ':', '"', '/', '\\', '|', '?', '*', '#', '%':
		return true
	}
	return r <= 0x1f || r == 0x7f
}

// SafeFmt 将名称中的不安全字符、"." 与 "-" 替换为 "_"，折叠连续下划线，
// 并去除首尾空白与下划线；结果为空时返回 fallback。
func SafeFmt(name, fallback string) string {
	var b strings.Builder
	b.Grow(len(name))
	for _, r := range name {
		if osUnsafe(r) {
			b.WriteByte('_')
			continue
		}
		switch r {
		case '.', '-':
			b.WriteByte('_')
		default:
			b.WriteRune(r)
		}
	}
	out := b.String()
	for strings.Contains(out, "__") {
		out = strings.ReplaceAll(out, "__", "_")
	}
	out = strings.Trim(strings.TrimSpace(out), "_")
	if out == "" {
		return fallback
	}
	return out
}

// SafeTruncate 在 UTF-8 字节长度超过 maxLen 时，保留前 3 个字符与尾部字符，
// 中间以 "..." 连接，确保结果不超过 maxLen 字节；未超限时原样返回。
// 第二个返回值表示是否发生了截断。
func SafeTruncate(name string, maxLen int) (string, bool, error) {
	if len(name) <= maxLen {
		return name, false, nil
	}
	if maxLen < minTruncateLen {
		return "", false, fmt.Errorf(
			"max_length (%d) is too small for SafeTruncate, must be >= %d", maxLen, minTruncateLen)
	}

	runes := []rune(name)
	head := runes
	if len(head) > truncHead {
		head = head[:truncHead]
	}
	// 按码点切片，长度按 UTF-8 字节数校验。
	tailLen := maxLen - truncHead - len(truncEllipsis)
	// 尾部长度超过字符数时取整个名称，等价于 Python 的 name[-tail_len:]。
	if tailLen > len(runes) {
		tailLen = len(runes)
	}
	trunc := string(head) + truncEllipsis + string(runes[len(runes)-tailLen:])
	for len(trunc) > maxLen && tailLen > 1 {
		tailLen--
		trunc = string(head) + truncEllipsis + string(runes[len(runes)-tailLen:])
	}
	if len(trunc) > maxLen {
		return "", false, fmt.Errorf(
			"max_length (%d) is too small for SafeTruncate: the first %d characters already exceed %d bytes",
			maxLen, truncHead, maxLen-truncHead)
	}
	return trunc, true, nil
}
