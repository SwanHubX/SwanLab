package fs

import (
	"path/filepath"
	"strings"
	"testing"
)

// 期望值与 Python SDK 的 fs.safe_fmt 一致。
func TestSafeFmt(t *testing.T) {
	tests := []struct {
		in, want string
	}{
		{"abc<>:\"/\\|?*#%def", "abc_def"},   // 连续非法字符折叠为单个下划线
		{"a.b-c_d", "a_b_c_d"},               // . 与 - 替换
		{"...---", "unknown"},                // 全非法 → fallback
		{"实验-run.01", "实验_run_01"},           // 多字节字符保留
		{"a//b", "a_b"},                      // 斜杠清洗
		{"  spaces  ", "spaces"},             // 首尾空白剥离
		{"__x__", "x"},                       // 折叠 + 剥离
		{"a\x00b\x7fc", "a_b_c"},             // 控制字符与 DEL
		{"normal_id-0123", "normal_id_0123"}, // 常规 run_id
	}
	for _, tc := range tests {
		if got := SafeFmt(tc.in, "unknown"); got != tc.want {
			t.Fatalf("SafeFmt(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

func TestSafeTruncate(t *testing.T) {
	tests := []struct {
		in      string
		max     int
		want    string
		trunced bool
	}{
		{"abcdefghijklmnopqrstuvwxyz", 10, "abc...wxyz", true},
		{"abcdefghijklmnopqrstuvwxyz", 7, "abc...z", true},
		{"abcdefghijklmnopqrstuvwxyz", 8, "abc...yz", true},
		{"short", 10, "short", false},
		{strings.Repeat("a", 99), 100, strings.Repeat("a", 99), false},
		// 101 字节触发截断：3 + "..." + 94 字节
		{strings.Repeat("a", 101), 100, "aaa" + "..." + strings.Repeat("a", 94), true},
		// 多字节：头部 3 字符 9 字节，尾部按码点收缩至总长 ≤ 100
		{strings.Repeat("实验", 80), 100, "实验实" + "..." + strings.Repeat("验实", 14) + "验", true},
		// 极端多字节：é（2 字节）截断到 20
		{strings.Repeat("é", 60), 20, "ééé" + "..." + strings.Repeat("é", 5), true},
	}
	for _, tc := range tests {
		got, trunced, err := SafeTruncate(tc.in, tc.max)
		if err != nil {
			t.Fatalf("SafeTruncate(%q, %d) 报错: %v", tc.in, tc.max, err)
		}
		if got != tc.want || trunced != tc.trunced {
			t.Fatalf("SafeTruncate(%q, %d) = (%q, %v), want (%q, %v)", tc.in, tc.max, got, trunced, tc.want, tc.trunced)
		}
		if len(got) > tc.max {
			t.Fatalf("结果 %d 字节超过上限 %d", len(got), tc.max)
		}
	}

	// max < 7 且确实需要截断时报错；不需要截断时不报错。
	if _, _, err := SafeTruncate("abcdefgh", 3); err == nil {
		t.Fatal("max=3 截断应报错")
	}
	if got, _, err := SafeTruncate("abc", 3); err != nil || got != "abc" {
		t.Fatalf("短名不应报错: %q %v", got, err)
	}

	// 头部不足 3 个字符时不能 panic，仍按超限返回错误。
	if _, _, err := SafeTruncate("🚀🚀", 7); err == nil {
		t.Fatal("多字节短名截断应报错")
	}

	// maxLen - 6 超过字符数时，尾部取整个名称，不能因负下标 panic。
	// 26 个 4 字节 emoji 共 104 字节，maxLen=100 触发截断且结果不超过上限。
	got, trunced, err := SafeTruncate(strings.Repeat("🚀", 26), 100)
	if err != nil || !trunced {
		t.Fatalf("多字节长名 = (%q, %v, %v)", got, trunced, err)
	}
	if n := len(got); n > 100 {
		t.Fatalf("多字节长名截断后 %d 字节超过上限", n)
	}
}

// TestWindowsParts 覆盖 Windows 语义的路径拆分。
func TestWindowsParts(t *testing.T) {
	cases := []struct {
		in   string
		want []string
	}{
		{`sub\model.py`, []string{"sub", "model.py"}},
		{`a/b\c`, []string{"a", "b", "c"}},
		{`./a//b`, []string{"a", "b"}},
		{`..\a`, []string{"..", "a"}},
		{"", nil},
	}
	for _, tc := range cases {
		got := WindowsParts(tc.in)
		if len(got) != len(tc.want) {
			t.Fatalf("WindowsParts(%q) = %v, want %v", tc.in, got, tc.want)
		}
		for i := range got {
			if got[i] != tc.want[i] {
				t.Fatalf("WindowsParts(%q) = %v, want %v", tc.in, got, tc.want)
			}
		}
	}
}

// TestMirrorPath 覆盖镜像路径与 Sender 回退路径共用的分隔符语义。
func TestMirrorPath(t *testing.T) {
	base := filepath.Join("run", "files")
	cases := []struct {
		name string
		want string
	}{
		{`sub\model.py`, filepath.Join(base, "sub", "model.py")},
		{`a/b\c`, filepath.Join(base, "a", "b", "c")},
		{`./a//b`, filepath.Join(base, "a", "b")},
	}
	for _, tc := range cases {
		if got := MirrorPath(base, tc.name); got != tc.want {
			t.Fatalf("MirrorPath(%q) = %q, want %q", tc.name, got, tc.want)
		}
	}
}

// TestWindowsBase 覆盖 Windows 语义的文件名提取。
func TestWindowsBase(t *testing.T) {
	cases := []struct {
		in   string
		want string
	}{
		{`C:\runs\files\config.yaml`, "config.yaml"},
		{"a/b", "b"},
		{"plain", "plain"},
		{"", ""},
	}
	for _, tc := range cases {
		if got := WindowsBase(tc.in); got != tc.want {
			t.Fatalf("WindowsBase(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}
