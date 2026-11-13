package mime

import "testing"

// TestTypeMatchesPythonRegistry 覆盖与 Python 自定义注册表逐项对齐的扩展名。
func TestTypeMatchesPythonRegistry(t *testing.T) {
	cases := []struct {
		path string
		want string
	}{
		{"config.yaml", "application/yaml"},
		{"config.yml", "application/yaml"},
		{"photo.png", "image/png"},
		{"photo.jpg", "image/jpeg"},
		{"photo.jpeg", "image/jpeg"},
		{"anim.gif", "image/gif"},
		{"photo.bmp", "image/bmp"},
		{"photo.webp", "image/webp"},
		{"logo.svg", "image/svg+xml"},
		{"logo.svgz", "image/svg+xml"},
		{"scan.tiff", "image/tiff"},
		{"scan.tif", "image/tiff"},
		{"train.py", "text/x-python"},
		{"notes.txt", "text/plain"},
		{"meta.json", "application/json"},
		{"README.md", "text/markdown"},
		{"README.markdown", "text/markdown"},
	}
	for _, tc := range cases {
		if got := Type(tc.path); got != tc.want {
			t.Errorf("Type(%q) = %q, want %q", tc.path, got, tc.want)
		}
	}
}

// TestTypeIgnoresCase 覆盖扩展名大小写不敏感。
func TestTypeIgnoresCase(t *testing.T) {
	for path, want := range map[string]string{
		"PHOTO.PNG":   "image/png",
		"Config.YAML": "application/yaml",
		"NOTES.TXT":   "text/plain",
	} {
		if got := Type(path); got != want {
			t.Errorf("Type(%q) = %q, want %q", path, got, want)
		}
	}
}

// TestTypeFallsBack 覆盖未知扩展名、无扩展名与含目录的路径。
func TestTypeFallsBack(t *testing.T) {
	cases := []struct {
		path string
		want string
	}{
		{"weights.swanlab-unknown", DefaultType},
		{"README", DefaultType},
		{"", DefaultType},
		{"/tmp/dir.v2/photo.JPG", "image/jpeg"}, // 目录里的点不能被当成扩展名
		{"C:\\runs\\exp\\train.py", "text/x-python"},
	}
	for _, tc := range cases {
		if got := Type(tc.path); got != tc.want {
			t.Errorf("Type(%q) = %q, want %q", tc.path, got, tc.want)
		}
	}
}
