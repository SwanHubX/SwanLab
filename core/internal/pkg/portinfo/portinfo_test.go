package portinfo

import (
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
)

func validUnixContent() string {
	return "protocol=1\npid=4242\nunix=/tmp/swanlab/core.sock\nEOF\n"
}

func TestMarshalUnixFormat(t *testing.T) {
	data, err := Marshal(&Info{Protocol: 1, PID: 4242, UnixPath: "/tmp/swanlab/core.sock"})
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if got := string(data); got != validUnixContent() {
		t.Fatalf("Marshal output mismatch:\n got: %q\nwant: %q", got, validUnixContent())
	}
}

func TestMarshalSockFormat(t *testing.T) {
	data, err := Marshal(&Info{Protocol: 1, PID: 4242, SockPort: 12345})
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	want := "protocol=1\npid=4242\nsock=12345\nEOF\n"
	if got := string(data); got != want {
		t.Fatalf("Marshal output mismatch:\n got: %q\nwant: %q", got, want)
	}
}

func TestWriteAndParseRoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "core.port")
	want := Info{Protocol: 1, PID: 4242, UnixPath: "/tmp/swanlab/core.sock"}
	err := WriteFile(path, &want)
	if err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	got, err := ParseFile(path)
	if err != nil {
		t.Fatalf("ParseFile: %v", err)
	}
	if got != want {
		t.Fatalf("round trip mismatch: got %+v, want %+v", got, want)
	}
	if runtime.GOOS != "windows" {
		var info os.FileInfo
		info, err = os.Stat(path)
		if err != nil {
			t.Fatalf("Stat: %v", err)
		}
		if perm := info.Mode().Perm(); perm != filePerm {
			t.Fatalf("port-file perm = %o, want %o", perm, filePerm)
		}
	}
	// 原子写入不应遗留临时文件
	var entries []os.DirEntry
	entries, err = os.ReadDir(dir)
	if err != nil {
		t.Fatalf("ReadDir: %v", err)
	}
	if len(entries) != 1 {
		t.Fatalf("runtime dir has %d entries after write, want 1", len(entries))
	}
}

func TestWriteFileAtomicallyReplaces(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "core.port")
	first := Info{Protocol: 1, PID: 4242, UnixPath: "/tmp/a.sock"}
	if err := WriteFile(path, &first); err != nil {
		t.Fatalf("first WriteFile: %v", err)
	}
	second := Info{Protocol: 1, PID: 4343, UnixPath: "/tmp/b.sock"}
	if err := WriteFile(path, &second); err != nil {
		t.Fatalf("second WriteFile: %v", err)
	}
	got, err := ParseFile(path)
	if err != nil {
		t.Fatalf("ParseFile: %v", err)
	}
	if got != second {
		t.Fatalf("after replace got %+v, want %+v", got, second)
	}
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("ReadDir: %v", err)
	}
	if len(entries) != 1 {
		t.Fatalf("runtime dir has %d entries after replace, want 1", len(entries))
	}
}

func TestWriteFileRejectsInvalidInfo(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "core.port")
	cases := map[string]Info{
		"bad protocol":  {Protocol: 2, PID: 4242, UnixPath: "/tmp/a.sock"},
		"both endings":  {Protocol: 1, PID: 4242, UnixPath: "/tmp/a.sock", SockPort: 80},
		"no endpoint":   {Protocol: 1, PID: 4242},
		"relative path": {Protocol: 1, PID: 4242, UnixPath: "tmp/a.sock"},
		"port range":    {Protocol: 1, PID: 4242, SockPort: 65536},
		"pid zero":      {Protocol: 1, PID: 0, UnixPath: "/tmp/a.sock"},
		"pid negative":  {Protocol: 1, PID: -1, UnixPath: "/tmp/a.sock"},
		"pid over max":  {Protocol: 1, PID: maxPID + 1, UnixPath: "/tmp/a.sock"},
	}
	for name, info := range cases {
		if err := WriteFile(path, &info); err == nil {
			t.Fatalf("%s: WriteFile unexpectedly succeeded", name)
		}
	}
	if _, err := os.Stat(path); !os.IsNotExist(err) {
		t.Fatalf("invalid info must not create target file, stat err = %v", err)
	}
	entries, err := os.ReadDir(dir)
	if err != nil {
		t.Fatalf("ReadDir: %v", err)
	}
	if len(entries) != 0 {
		t.Fatalf("invalid info left %d temp files, want 0", len(entries))
	}
}

func TestParseAcceptsValidVariants(t *testing.T) {
	// EOF 行末尾可无换行
	content := strings.TrimSuffix(validUnixContent(), "\n")
	if _, err := Parse([]byte(content)); err != nil {
		t.Fatalf("parse without trailing newline: %v", err)
	}
	// 键值行顺序不影响解析
	reordered := "unix=/tmp/swanlab/core.sock\npid=4242\nprotocol=1\nEOF\n"
	if _, err := Parse([]byte(reordered)); err != nil {
		t.Fatalf("parse reordered lines: %v", err)
	}
	// sock 端点（Windows 或 UDS 回退形态）
	sockContent := "protocol=1\npid=4242\nsock=12345\nEOF\n"
	info, err := Parse([]byte(sockContent))
	if err != nil {
		t.Fatalf("parse sock content: %v", err)
	}
	if info.SockPort != 12345 || info.UnixPath != "" {
		t.Fatalf("sock parse got %+v", info)
	}
}

func TestParseRejectsMalformed(t *testing.T) {
	cases := map[string]string{
		"missing EOF":          "protocol=1\npid=4242\nunix=/tmp/a.sock\n",
		"EOF not own line":     "protocol=1\npid=4242\nunix=/tmp/a.sockEOF\n",
		"content after EOF":    validUnixContent() + "extra\n",
		"torn write":           "protocol=1\npid=4242\n",
		"empty body":           "EOF\n",
		"empty line":           "protocol=1\n\npid=4242\nunix=/tmp/a.sock\nEOF\n",
		"not key value":        "protocol=1\npid\nunix=/tmp/a.sock\nEOF\n",
		"unknown key":          "protocol=1\npid=4242\nunix=/tmp/a.sock\nextra=1\nEOF\n",
		"legacy mode key":      "protocol=1\npid=4242\nmode=owner\nunix=/tmp/a.sock\nEOF\n",
		"legacy auth key":      "protocol=1\npid=4242\nunix=/tmp/a.sock\nauth_token=abc\nEOF\n",
		"duplicate key":        "protocol=1\nprotocol=1\npid=4242\nunix=/tmp/a.sock\nEOF\n",
		"unknown protocol":     "protocol=2\npid=4242\nunix=/tmp/a.sock\nEOF\n",
		"non numeric protocol": "protocol=abc\npid=4242\nunix=/tmp/a.sock\nEOF\n",
		"missing protocol":     "pid=4242\nunix=/tmp/a.sock\nEOF\n",
		"missing pid":          "protocol=1\nunix=/tmp/a.sock\nEOF\n",
		"missing endpoint":     "protocol=1\npid=4242\nEOF\n",
		"both endpoints":       "protocol=1\npid=4242\nunix=/tmp/a.sock\nsock=12345\nEOF\n",
		"relative unix path":   "protocol=1\npid=4242\nunix=tmp/a.sock\nEOF\n",
		"oversize unix path":   "protocol=1\npid=4242\nunix=/" + strings.Repeat("a", maxUnixPathLen) + "\nEOF\n",
		"port zero":            "protocol=1\npid=4242\nsock=0\nEOF\n",
		"port range":           "protocol=1\npid=4242\nsock=65536\nEOF\n",
		"port leading zero":    "protocol=1\npid=4242\nsock=01234\nEOF\n",
		"port not digits":      "protocol=1\npid=4242\nsock=12a45\nEOF\n",
		"pid zero":             "protocol=1\npid=0\nunix=/tmp/a.sock\nEOF\n",
		"pid negative":         "protocol=1\npid=-1\nunix=/tmp/a.sock\nEOF\n",
		"pid leading zero":     "protocol=1\npid=04242\nunix=/tmp/a.sock\nEOF\n",
		"pid not digits":       "protocol=1\npid=42a42\nunix=/tmp/a.sock\nEOF\n",
		"pid out of range":     "protocol=1\npid=2147483648\nunix=/tmp/a.sock\nEOF\n",
	}
	for name, content := range cases {
		if _, err := Parse([]byte(content)); err == nil {
			t.Fatalf("%s: Parse unexpectedly succeeded", name)
		}
	}
}

func TestParseRejectsOversizeContent(t *testing.T) {
	content := "protocol=1\npid=4242\nunix=/tmp/" + strings.Repeat("a", MaxFileSize) + "\nEOF\n"
	if _, err := Parse([]byte(content)); err == nil {
		t.Fatal("Parse unexpectedly succeeded for oversize content")
	}
}
