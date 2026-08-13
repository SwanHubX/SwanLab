package console

import (
	"io"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

func TestEmitPlainFormat(t *testing.T) {
	plain := emit("debug", ansiGrey, "hello world", "hello world", false)
	// 2026-03-14 10:23:45.124 | DEBUG    | console.TestEmitPlainFormat:line - hello world
	re := regexp.MustCompile(`^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3} \| DEBUG    \| .* - hello world$`)
	if !re.MatchString(plain) {
		t.Fatalf("unexpected plain format: %q", plain)
	}
}

func TestEmitLevelPadding(t *testing.T) {
	for _, lvl := range []struct{ name, padded string }{
		{"info", "INFO    "},
		{"warning", "WARNING "},
		{"error", "ERROR   "},
	} {
		plain := emit(lvl.name, "", "x", "x", false)
		if !strings.Contains(plain, "| "+lvl.padded+" |") {
			t.Fatalf("level %q not padded to %q in %q", lvl.name, lvl.padded, plain)
		}
	}
}

func TestDebugGatedByEnv(t *testing.T) {
	Reset()
	defer Reset()
	dir := t.TempDir()
	if err := Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	old := debugEnabled
	debugEnabled = false
	Debug("should not appear")
	debugEnabled = old
	data, _ := os.ReadFile(filepath.Join(dir, "debug-core.log"))
	if strings.Contains(string(data), "should not appear") {
		t.Fatalf("Debug wrote while disabled: %q", data)
	}
	debugEnabled = true
	Debug("should appear")
	data, _ = os.ReadFile(filepath.Join(dir, "debug-core.log"))
	if !strings.Contains(string(data), "should appear") {
		t.Fatalf("Debug did not write while enabled: %q", data)
	}
}

func TestInfoPrintsSwanlabPrefix(t *testing.T) {
	oldColor := colorEnabled
	colorEnabled = false
	defer func() { colorEnabled = oldColor }()
	out := captureStderr(func() { Info("hello", "world") })
	trimmed := strings.TrimSpace(out)
	if !strings.HasPrefix(trimmed, "swanlab:") {
		t.Fatalf("missing 'swanlab:' prefix in %q", out)
	}
	if !strings.Contains(trimmed, "hello world") {
		t.Fatalf("missing message in %q", out)
	}
}

func TestTraceNilIsNoOp(t *testing.T) {
	Reset()
	defer Reset()
	dir := t.TempDir()
	if err := Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	Trace(nil, "prefix")
	data, _ := os.ReadFile(filepath.Join(dir, "debug-core.log"))
	if strings.Contains(string(data), "prefix") {
		t.Fatalf("Trace(nil) wrote unexpectedly: %q", data)
	}
}

func TestTraceWritesErrorAndStack(t *testing.T) {
	Reset()
	defer Reset()
	dir := t.TempDir()
	if err := Init(dir); err != nil {
		t.Fatalf("Init: %v", err)
	}
	Trace(&testErr{"boom"}, "ctx")
	data, _ := os.ReadFile(filepath.Join(dir, "debug-core.log"))
	body := string(data)
	if !strings.Contains(body, "ctx: boom") {
		t.Fatalf("missing 'ctx: boom' in file: %q", body)
	}
	if !strings.Contains(body, "goroutine") {
		t.Fatalf("missing stack trace in file: %q", body)
	}
}

func TestTruthy(t *testing.T) {
	for _, v := range []string{"true", "TRUE", "1", "yes", "on", " On "} {
		if !truthy(v) {
			t.Fatalf("truthy(%q) = false, want true", v)
		}
	}
	for _, v := range []string{"", "false", "0", "no", "maybe"} {
		if truthy(v) {
			t.Fatalf("truthy(%q) = true, want false", v)
		}
	}
}

type testErr struct{ msg string }

func (e *testErr) Error() string { return e.msg }

// captureStderr 在 f 执行期间将 os.Stderr 重定向到管道以捕获输出。
func captureStderr(f func()) string {
	r, w, _ := os.Pipe()
	old := os.Stderr
	os.Stderr = w
	done := make(chan string)
	go func() { b, _ := io.ReadAll(r); done <- string(b) }()
	f()
	w.Close()
	os.Stderr = old
	return <-done
}
