// Package console 是 SwanLab core 面向用户的终端输出模块，并桥接到诊断日志子包 log。
// 镜像 swanlab/sdk/internal/pkg/console：
//   - info:    "swanlab: <message>"，蓝色前缀
//   - debug:   loguru 风格，灰色；仅 SWANLAB_DEBUG=true 时输出
//   - warning: loguru 风格，黄色
//   - error:   loguru 风格，红色
//   - trace:   error 变体，额外记录 goroutine 栈
//
// SWANLAB_DEBUG=true 时，所有级别都会同步镜像到诊断日志文件（经 log 子包）。
package console

import (
	"fmt"
	"os"
	"runtime"
	"runtime/debug"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console/log"
)

const name = "swanlab"

// ANSI 样式码（仅 colorEnabled 时生效）
const (
	ansiReset  = "\033[0m"
	ansiBold   = "\033[1m"
	ansiDim    = "\033[2m"
	ansiGrey   = "\033[38;5;242m" // rich "grey54" ≈ 256 色板 242
	ansiYellow = "\033[33m"
	ansiRed    = "\033[31m"
	ansiBlue   = "\033[34m"
)

var (
	mu           sync.Mutex // 串行化终端写入，保证日志行完整不被打断
	debugEnabled = truthy(os.Getenv("SWANLAB_DEBUG"))
	colorEnabled = detectColor()
	thisFile     string
)

func init() {
	// 记录本源码文件路径，供 callerLocation 跳过属于本文件的栈帧
	if _, file, _, ok := runtime.Caller(0); ok {
		thisFile = file
	}
}

// Init 绑定诊断日志目录（flush 内存缓冲）；dir 为空时禁用文件持久化。对应 Python console.init。
func Init(dir string) error { return log.Init(dir) }

// Reset 将诊断日志重置为初始的内存缓冲状态。
func Reset() { log.Reset() }

// Debug 输出 loguru 风格的调试日志。SWANLAB_DEBUG=false 时无操作。
//
//	2026-03-14 10:23:45.124 | DEBUG    | pkg.func:line - message
func Debug(args ...any) {
	if !debugEnabled {
		return
	}
	msg := joinArgs(args)
	log.Debug(emit("debug", ansiGrey, msg, msg, true))
}

// Debugf 是 Debug 的格式化变体。
func Debugf(format string, args ...any) { Debug(fmt.Sprintf(format, args...)) }

// Info 输出形如 "swanlab: <message>" 的通知。
func Info(args ...any) {
	msg := joinArgs(args)
	writeTTY(paint(name, ansiBlue+ansiBold) + ":" + " " + msg)
	if debugEnabled {
		log.Info(emit("info", "", msg, msg, false))
	}
}

// Infof 是 Info 的格式化变体。
func Infof(format string, args ...any) { Info(fmt.Sprintf(format, args...)) }

// Warning 输出 loguru 风格的警告（黄色）。
func Warning(args ...any) {
	msg := joinArgs(args)
	log.Warning(emit("warning", ansiYellow, msg, msg, true))
}

// Warningf 是 Warning 的格式化变体。
func Warningf(format string, args ...any) { Warning(fmt.Sprintf(format, args...)) }

// Error 输出 loguru 风格的错误（红色）。
func Error(args ...any) {
	msg := joinArgs(args)
	log.Error(emit("error", ansiRed, msg, msg, true))
}

// Errorf 是 Error 的格式化变体。
func Errorf(format string, args ...any) { Error(fmt.Sprintf(format, args...)) }

// Trace 记录 err 及其 goroutine 栈。终端只显示短消息，诊断日志文件额外记录完整栈。
// err 为 nil 时无操作。这是 Python console.trace（读取 sys.exc_info）的 Go 适配版。
func Trace(err error, args ...any) {
	if err == nil {
		return
	}
	prefix := joinArgs(args)
	stack := strings.TrimSpace(string(debug.Stack()))
	// 与 Python 回溯不同，debug.Stack() 不含错误信息，因此此处显式前置，
	// 保证文件记录自包含（终端仅显示短消息）
	shortMsg := err.Error()
	fileMsg := err.Error() + "\n" + stack
	if prefix != "" {
		shortMsg = prefix + ": " + shortMsg
		fileMsg = prefix + ": " + fileMsg
	}
	log.Error(emit("error", ansiRed, shortMsg, fileMsg, true))
}

// emit 构建并（可选）打印 loguru 风格日志行，返回纯文本（供诊断日志文件使用）。
// consoleMsg 用于终端显示，fileMsg 用于文件持久化（仅 Trace 时二者不同）。
func emit(level, color, consoleMsg, fileMsg string, writeToTTY bool) string {
	ts := now()
	loc := callerLocation()
	lvl := fmt.Sprintf("%-8s", strings.ToUpper(level))
	plain := ts + " | " + lvl + " | " + loc + " - " + fileMsg
	if !writeToTTY {
		return plain
	}
	if colorEnabled {
		sep := color + ansiDim
		writeTTY(paint(ts, color) + paint(" | ", sep) + paint(lvl, color+ansiBold) +
			paint(" | ", sep) + paint(loc, color) + paint(" - ", sep) + consoleMsg)
	} else {
		writeTTY(ts + " | " + lvl + " | " + loc + " - " + consoleMsg)
	}
	return plain
}

func now() string { return time.Now().Format("2006-01-02 15:04:05.000") }

// callerLocation 返回首个不属于本文件的栈帧的 "pkg.func:line"，对应 Python _caller_location 的栈遍历。
func callerLocation() string {
	pcs := make([]uintptr, 32)
	n := runtime.Callers(2, pcs)
	frames := runtime.CallersFrames(pcs[:n])
	for {
		fr, more := frames.Next()
		if fr.File != "" && fr.File != thisFile {
			return shortFunc(fr.Function) + ":" + strconv.Itoa(fr.Line)
		}
		if !more {
			break
		}
	}
	return name
}

func shortFunc(full string) string {
	if i := strings.LastIndex(full, "/"); i >= 0 {
		full = full[i+1:]
	}
	return full
}

// joinArgs 以单空格拼接 args，对应 rich print / _to_plain_text 的行为。
func joinArgs(args []any) string {
	parts := make([]string, len(args))
	for i, a := range args {
		parts[i] = fmt.Sprint(a)
	}
	return strings.Join(parts, " ")
}

// paint 在启用颜色时用 ANSI 码包裹 s，否则原样返回。
func paint(s, codes string) string {
	if !colorEnabled || codes == "" {
		return s
	}
	return codes + s + ansiReset
}

func writeTTY(s string) {
	mu.Lock()
	defer mu.Unlock()
	fmt.Fprintln(os.Stderr, s)
}

// detectColor 在 stderr 为终端且未设置 NO_COLOR 时启用 ANSI 颜色（https://no-color.org）。
// 使用标准库的字符设备检测，避免引入额外依赖。
func detectColor() bool {
	if os.Getenv("NO_COLOR") != "" {
		return false
	}
	if fi, err := os.Stderr.Stat(); err == nil && (fi.Mode()&os.ModeCharDevice) != 0 {
		return true
	}
	return false
}

// truthy 判断 v 是否为 SWANLAB_DEBUG 的真值。
func truthy(v string) bool {
	switch strings.ToLower(strings.TrimSpace(v)) {
	case "true", "1", "yes", "on":
		return true
	}
	return false
}
