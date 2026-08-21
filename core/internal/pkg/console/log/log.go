// Package log 是 SwanLab core 的内部诊断日志模块，镜像 swanlab/sdk/internal/pkg/console/log。
//
// 零额外依赖：日志在内存中缓冲，直到 Init 绑定日志目录后一次性 flush 到轮转文件
// (debug-core.log)，文件权限强制为 0600。
//
// 与 console 模块互补：
//   - console：面向用户的终端美化输出
//   - log：面向开发者/运维的文件持久化诊断日志
//
// console 模块负责按 loguru 风格预格式化每一行（时间戳 | 级别 | 位置 - 消息），
// 将纯文本字符串传入本模块；本模块原样写入，对应 Python logging.Formatter(fmt="%(message)s")。
//
// 采用 Go 标准库 log 包的混合模式：Logger 结构体承载状态，包级默认实例 std 提供便捷 API
// （Init/Debug 等转发给 std）。测试可直接 New() 构造独立实例，互不干扰、无需 Reset。
package log

import (
	"fmt"
	"os"
	"path/filepath"
	"sync"
)

const (
	memoryCapacity = 1024 // Init 绑定文件前的最大内存缓冲条数
	filename       = "debug-core.log"
)

// Logger 是一个可绑定到文件、支持轮转的诊断日志器。使用 New() 构造。
type Logger struct {
	mu      sync.Mutex
	dir     string   // Init 绑定的目录，绑定前为空
	file    *os.File // 当前 debug-core.log 句柄
	size    int64    // 当前文件大小（字节）
	buf     []record // 绑定前的内存缓冲
	bound   bool     // Init 是否已调用
	enabled bool     // 文件输出是否生效（Init 传入非空目录）

	maxBytes    int64 // 单文件上限，超过则轮转（0 表示不限制）
	backupCount int   // 保留的备份数
}

// record 是等待绑定文件的缓冲日志条目
type record struct{ msg string }

// New 返回一个未绑定文件的 Logger（内存缓冲态），使用默认轮转配置（10 MB × 3 备份）。
// 按单次实验生命周期设计，不重复挂载
func New() *Logger {
	return &Logger{
		maxBytes:    10 * 1024 * 1024, // 单文件上限 10 MB
		backupCount: 3,                // 保留 debug-core.log.1 .. .3（共 ~40 MB）
	}
}

// Init 绑定诊断日志目录（将内存缓冲 flush 到 debug-core.log）；dir 为空时禁用持久化并丢弃缓冲。
// 重复调用为空操作，需先 Reset 才能重新绑定。
func (l *Logger) Init(dir string) error {
	l.mu.Lock()
	defer l.mu.Unlock()
	if l.bound {
		return nil
	}
	if dir == "" {
		l.buf = nil
		l.enabled = false
		l.bound = true
		return nil
	}
	if info, err := os.Stat(dir); err != nil || !info.IsDir() {
		return fmt.Errorf("log directory does not exist: %s", dir)
	}
	if err := l.openLocked(dir); err != nil {
		return err
	}
	for _, r := range l.buf {
		l.writeLocked(r)
	}
	l.buf = nil
	return nil
}

// openLocked 在 dir 下打开 debug-core.log 并记录当前文件大小。
// 调用方须持有 mu。POSIX 下强制 0600 权限（Windows 下尽力而为）。
func (l *Logger) openLocked(dir string) error {
	path := filepath.Join(dir, filename)
	f, err := os.OpenFile(path, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0o600)
	if err != nil {
		return err
	}
	_ = os.Chmod(path, 0o600) // 强制 0600，不受 umask 影响
	info, err := f.Stat()
	if err != nil {
		f.Close()
		return err
	}
	l.dir = dir
	l.file = f
	l.size = info.Size()
	l.enabled = true
	l.bound = true
	return nil
}

// Reset 将日志器重置为初始的内存缓冲状态，关闭已打开的文件句柄。
// 这里有个设计细节，因为 Init 绑定后 Logger 对象获得了文件句柄，但是 Go 没有析构函数，所以Logger被释放前需要手动调用Reset
// 不过这并不是什么大问题，因为core在设计上对应一次实验生命周期，所以FD泄漏的风险不大，且Reset在core的生命周期结束时会被调用。
func (l *Logger) Reset() {
	l.mu.Lock()
	defer l.mu.Unlock()
	if l.file != nil {
		l.file.Close()
		l.file = nil
	}
	l.dir = ""
	l.size = 0
	l.buf = nil
	l.enabled = false
	l.bound = false
}

// Debug/Info/Warning/Error/Critical 均原样写入 msg（级别名义化，
// 对应 Python 仅输出消息的 formatter）；保留不同方法名以与 Python log 模块 API 对齐。
func (l *Logger) Debug(msg string)    { l.logLine(msg) }
func (l *Logger) Info(msg string)     { l.logLine(msg) }
func (l *Logger) Warning(msg string)  { l.logLine(msg) }
func (l *Logger) Error(msg string)    { l.logLine(msg) }
func (l *Logger) Critical(msg string) { l.logLine(msg) }

func (l *Logger) logLine(msg string) {
	l.mu.Lock()
	defer l.mu.Unlock()
	if !l.bound {
		// 阶段 1：缓冲到内存，直到 Init 绑定文件（超出容量时丢弃最旧记录）
		l.buf = append(l.buf, record{msg: msg})
		if len(l.buf) > memoryCapacity {
			l.buf = l.buf[len(l.buf)-memoryCapacity:]
		}
		return
	}
	if !l.enabled {
		return // 禁用模式（Init("")）：丢弃
	}
	l.writeLocked(record{msg: msg})
}

// writeLocked 将一条记录追加到 debug-core.log，若超过 maxBytes 则先轮转。调用方须持有 mu。
func (l *Logger) writeLocked(r record) {
	if l.file == nil {
		return
	}
	line := r.msg + "\n"
	if l.maxBytes > 0 && l.size+int64(len(line)) > l.maxBytes {
		l.rotateLocked()
	}
	if l.file == nil {
		return
	}
	if n, err := l.file.WriteString(line); err == nil {
		l.size += int64(n)
	}
}

// rotateLocked 移动备份（debug-core.log.1->.2->.3，丢弃最旧的），并打开新的 debug-core.log。调用方须持有 mu。
func (l *Logger) rotateLocked() {
	if l.file != nil {
		l.file.Close()
		l.file = nil
	}
	base := filepath.Join(l.dir, filename)
	for i := l.backupCount - 1; i >= 1; i-- {
		src := fmt.Sprintf("%s.%d", base, i)
		dst := fmt.Sprintf("%s.%d", base, i+1)
		_ = os.Remove(dst)
		_ = os.Rename(src, dst)
	}
	_ = os.Rename(base, base+".1")
	f, err := os.OpenFile(base, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0o600)
	if err != nil {
		return
	}
	_ = os.Chmod(base, 0o600)
	l.file = f
	l.size = 0
}

// -----------------------------------------------------------------------------
// 包级默认实例：便捷 API 转发给 std，对应标准库 log 包的设计。
// console 等只需要进程级单例的调用方直接用这些函数即可；需要隔离的测试用 New()。
// -----------------------------------------------------------------------------

var std = New()

// Init 在默认实例上调用 Init。详见 (*Logger).Init。
func Init(dir string) error { return std.Init(dir) }

// Reset 在默认实例上调用 Reset。详见 (*Logger).Reset。
func Reset() { std.Reset() }

// Debug 在默认实例上调用 Debug。详见 (*Logger).Debug。
func Debug(msg string) { std.Debug(msg) }

// Info 在默认实例上调用 Info。详见 (*Logger).Info。
func Info(msg string) { std.Info(msg) }

// Warning 在默认实例上调用 Warning。详见 (*Logger).Warning。
func Warning(msg string) { std.Warning(msg) }

// Error 在默认实例上调用 Error。详见 (*Logger).Error。
func Error(msg string) { std.Error(msg) }

// Critical 在默认实例上调用 Critical。详见 (*Logger).Critical。
func Critical(msg string) { std.Critical(msg) }
