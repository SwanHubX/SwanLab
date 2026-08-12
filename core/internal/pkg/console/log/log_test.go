package log

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// tempLogDir 返回一个临时目录，并注册清理：在 t.TempDir 的 RemoveAll 之前关闭 Logger 句柄。
//
// Logger 设计上不要求调用方调用 Reset——core 生命周期内会统一管理，FD 泄漏风险很低。
// 此处仅为测试场景：Windows 不允许删除仍被打开的文件，若不关闭句柄，t.TempDir 的
// RemoveAll 清理会失败（"being used by another process"）；POSIX 下删除已打开文件本就允许。
//
// 顺序要点：t.Cleanup 按注册逆序（LIFO）执行，而 t.TempDir 的 RemoveAll 也是一次 Cleanup。
// 这里先建目录、再注册 l.Reset，使 Reset 晚注册而先执行（先于 RemoveAll）。
func tempLogDir(t *testing.T, l *Logger) string {
	t.Helper()
	dir := t.TempDir()
	t.Cleanup(l.Reset)
	return dir
}

func TestBufferUntilInit(t *testing.T) {
	l := New()
	l.Debug("buffered line")
	dir := tempLogDir(t, l)
	if err := l.Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	data, err := os.ReadFile(filepath.Join(dir, filename))
	if err != nil {
		t.Fatalf("read file: %v", err)
	}
	if !strings.Contains(string(data), "buffered line") {
		t.Fatalf("buffered line not flushed to file; got %q", data)
	}
}

func TestInitDisabledDiscards(t *testing.T) {
	l := New()
	l.Debug("will be discarded")
	if err := l.Init(""); err != nil {
		t.Fatalf("Init: %v", err)
	}
	l.Debug("after disable") // 已禁用，应被丢弃且不 panic
}

func TestIdempotentInit(t *testing.T) {
	l := New()
	dir := tempLogDir(t, l)
	if err := l.Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	dir2 := t.TempDir()
	if err := l.Init(dir2); err != nil {
		t.Fatalf("second Init: %v", err)
	}
	// 第二次 Init 为空操作：dir2 不应包含日志文件
	if _, err := os.Stat(filepath.Join(dir2, filename)); err == nil {
		t.Fatalf("second Init should be a no-op, but file exists in dir2")
	}
}

func TestRotation(t *testing.T) {
	l := New()
	l.maxBytes = 50 // 调小上限以触发轮转
	dir := tempLogDir(t, l)
	if err := l.Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	for range 20 {
		l.Error("rotation-test-line")
	}
	if _, err := os.Stat(filepath.Join(dir, filename+".1")); err != nil {
		t.Fatalf("expected backup .1 after rotation: %v", err)
	}
}

func TestResetAllowsRebind(t *testing.T) {
	l := New()
	dir := t.TempDir()
	if err := l.Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	l.Reset()
	dir2 := tempLogDir(t, l)
	l.Debug("buffered after reset")
	if err := l.Init(dir2); err != nil {
		t.Fatalf("rebind Init: %v", err)
	}
	data, _ := os.ReadFile(filepath.Join(dir2, filename))
	if !strings.Contains(string(data), "buffered after reset") {
		t.Fatalf("expected rebind flush; got %q", data)
	}
}

func TestInitRejectsMissingDir(t *testing.T) {
	l := New()
	missing := filepath.Join(t.TempDir(), "does-not-exist")
	if err := l.Init(missing); err == nil {
		t.Fatalf("Init with missing dir should error")
	}
}

// TestPackageLevelWrappers 验证包级便捷函数确实转发给默认实例 std。
func TestPackageLevelWrappers(t *testing.T) {
	Reset()
	defer Reset()
	dir := t.TempDir()
	if err := Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	Debug("via-std")
	data, _ := os.ReadFile(filepath.Join(dir, filename))
	if !strings.Contains(string(data), "via-std") {
		t.Fatalf("package-level Debug did not forward to std; got %q", data)
	}
}
