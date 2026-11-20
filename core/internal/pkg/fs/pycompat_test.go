package fs

import (
	"strings"
	"testing"
)

// 期望值由 Python SDK 的 fs.safe_fmt / fs.safe_truncate 生成，用于锁定跨语言行为。
func TestSafeFmtMatchesPython(t *testing.T) {
	tests := []struct {
		in, fallback, want string
	}{
		{"", "unknown", "unknown"},
		{" ", "x", "x"},
		{"___", "", ""},
		{"a__b", "unknown", "a_b"},
		{"a.-b", "unknown", "a_b"},
		{".hidden", "unknown", "hidden"},
		{"a/b", "unknown", "a_b"},
		{"a\\b", "unknown", "a_b"},
		{"a:b", "unknown", "a_b"},
		{"a b", "unknown", "a b"},
		{" a ", "unknown", "a"},
		{"\ta\n", "unknown", "a"},
		{"\x00", "unknown", "unknown"},
		{"\x7f", "unknown", "unknown"},
		{"\u0085", "unknown", "unknown"},
		{"\u00a0", "unknown", "unknown"},
		{"\u3000", "unknown", "unknown"},
		{"a\u0301", "unknown", "a\u0301"},
		{"实验 run.01", "unknown", "实验 run_01"},
		{"🚀🚀", "unknown", "🚀🚀"},
		{"abc<>:\"/\\|?*#%def", "unknown", "abc_def"},
		{"...---", "unknown", "unknown"},
		{"  __ a __  ", "unknown", " a "},
		{strings.Repeat("a", 101), "unknown", strings.Repeat("a", 101)},
	}
	for _, tc := range tests {
		if got := SafeFmt(tc.in, tc.fallback); got != tc.want {
			t.Fatalf("SafeFmt(%q, %q) = %q, want %q", tc.in, tc.fallback, got, tc.want)
		}
	}
}

func TestSafeTruncateMatchesPython(t *testing.T) {
	tests := []struct {
		in      string
		maxLen  int
		want    string
		trunced bool
		wantErr bool
	}{
		{"abcdefghijklmnopqrstuvwxyz", 10, "abc...wxyz", true, false},
		{"abcdefghijklmnopqrstuvwxyz", 7, "abc...z", true, false},
		{"abcdefghijklmnopqrstuvwxyz", 8, "abc...yz", true, false},
		{"short", 10, "short", false, false},
		{"", 0, "", false, false},
		{"a", 0, "", false, true},
		{"a", 6, "a", false, false},
		{"ab", 6, "ab", false, false},
		{"实验", 5, "", false, true},
		{"实验", 6, "实验", false, false},
		{"node_modules", 8, "nod...es", true, false},
		{"中文指标名称", 20, "中文指标名称", false, false},
		{"abcdefg", 7, "abcdefg", false, false},
		{"abcdefgh", 7, "abc...h", true, false},
		{strings.Repeat("a", 99), 100, strings.Repeat("a", 99), false, false},
		{strings.Repeat("a", 101), 100, "aaa" + "..." + strings.Repeat("a", 94), true, false},
		{strings.Repeat("a", 200), 100, "aaa" + "..." + strings.Repeat("a", 94), true, false},
		{strings.Repeat("实验", 80), 100, "实验实" + "..." + strings.Repeat("验实", 14) + "验", true, false},
		{strings.Repeat("é", 60), 20, "ééé" + "..." + strings.Repeat("é", 5), true, false},
		{"🚀🚀", 7, "", false, true},
		{strings.Repeat("🚀", 26), 100, "🚀🚀🚀" + "..." + strings.Repeat("🚀", 21), true, false},
		{strings.Repeat("🚀", 26), 50, "🚀🚀🚀" + "..." + strings.Repeat("🚀", 8), true, false},
	}
	for _, tc := range tests {
		got, trunced, err := SafeTruncate(tc.in, tc.maxLen)
		if (err != nil) != tc.wantErr {
			t.Fatalf("SafeTruncate(%q, %d) err = %v, wantErr %v", tc.in, tc.maxLen, err, tc.wantErr)
		}
		if tc.wantErr {
			continue
		}
		if got != tc.want || trunced != tc.trunced {
			t.Fatalf("SafeTruncate(%q, %d) = (%q, %v), want (%q, %v)", tc.in, tc.maxLen, got, trunced, tc.want, tc.trunced)
		}
	}
}
