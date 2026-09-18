package portinfo

import (
	"bytes"
	"encoding/base64"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
)

// testToken 返回确定性的合法 256-bit token，避免测试依赖随机数。
func testToken() string {
	return base64.RawURLEncoding.EncodeToString(bytes.Repeat([]byte{0x42}, authTokenBytes))
}

func validUnixContent() string {
	return "protocol=1\nunix=/tmp/swanlab/core.sock\nauth=" + testToken() + "\nEOF\n"
}

func TestMarshalUnixFormat(t *testing.T) {
	data, err := Marshal(&Info{Protocol: 1, UnixPath: "/tmp/swanlab/core.sock", AuthToken: testToken()})
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if got := string(data); got != validUnixContent() {
		t.Fatalf("Marshal output mismatch:\n got: %q\nwant: %q", got, validUnixContent())
	}
}

func TestMarshalSockFormat(t *testing.T) {
	data, err := Marshal(&Info{Protocol: 1, SockPort: 12345, AuthToken: testToken()})
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	want := "protocol=1\nsock=12345\nauth=" + testToken() + "\nEOF\n"
	if got := string(data); got != want {
		t.Fatalf("Marshal output mismatch:\n got: %q\nwant: %q", got, want)
	}
}

func TestWriteAndParseRoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "core.port")
	token, err := NewAuthToken()
	if err != nil {
		t.Fatalf("NewAuthToken: %v", err)
	}
	want := Info{Protocol: 1, UnixPath: "/tmp/swanlab/core.sock", AuthToken: token}
	err = WriteFile(path, &want)
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
	first := Info{Protocol: 1, UnixPath: "/tmp/a.sock", AuthToken: testToken()}
	if err := WriteFile(path, &first); err != nil {
		t.Fatalf("first WriteFile: %v", err)
	}
	secondToken, err := NewAuthToken()
	if err != nil {
		t.Fatalf("NewAuthToken: %v", err)
	}
	second := Info{Protocol: 1, UnixPath: "/tmp/b.sock", AuthToken: secondToken}
	err = WriteFile(path, &second)
	if err != nil {
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
		"bad protocol":  {Protocol: 2, UnixPath: "/tmp/a.sock", AuthToken: testToken()},
		"both endings":  {Protocol: 1, UnixPath: "/tmp/a.sock", SockPort: 80, AuthToken: testToken()},
		"no endpoint":   {Protocol: 1, AuthToken: testToken()},
		"relative path": {Protocol: 1, UnixPath: "tmp/a.sock", AuthToken: testToken()},
		"port range":    {Protocol: 1, SockPort: 65536, AuthToken: testToken()},
		"bad token":     {Protocol: 1, UnixPath: "/tmp/a.sock", AuthToken: "short"},
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

func TestNewAuthToken(t *testing.T) {
	first, err := NewAuthToken()
	if err != nil {
		t.Fatalf("NewAuthToken: %v", err)
	}
	second, err := NewAuthToken()
	if err != nil {
		t.Fatalf("NewAuthToken: %v", err)
	}
	if first == second {
		t.Fatal("two generated tokens must differ")
	}
	if raw, err := base64.RawURLEncoding.Strict().DecodeString(first); err != nil || len(raw) != authTokenBytes {
		t.Fatalf("generated token is not 256-bit base64url: decode err = %v", err)
	}
}

func TestParseAcceptsValidVariants(t *testing.T) {
	// EOF 行末尾无换行同样合法
	content := strings.TrimSuffix(validUnixContent(), "\n")
	if _, err := Parse([]byte(content)); err != nil {
		t.Fatalf("parse without trailing newline: %v", err)
	}
	// 键值行顺序不影响解析
	reordered := "auth=" + testToken() + "\nunix=/tmp/swanlab/core.sock\nprotocol=1\nEOF\n"
	if _, err := Parse([]byte(reordered)); err != nil {
		t.Fatalf("parse reordered lines: %v", err)
	}
	// sock 端点（Windows 形态）
	sockContent := "protocol=1\nsock=12345\nauth=" + testToken() + "\nEOF\n"
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
		"missing EOF":          "protocol=1\nunix=/tmp/a.sock\nauth=" + testToken() + "\n",
		"EOF not own line":     "protocol=1\nunix=/tmp/a.sock\nauth=" + testToken() + "EOF\n",
		"content after EOF":    validUnixContent() + "extra\n",
		"torn write":           "protocol=1\nunix=/tmp/a.sock\n",
		"empty body":           "EOF\n",
		"empty line":           "protocol=1\n\nunix=/tmp/a.sock\nauth=" + testToken() + "\nEOF\n",
		"not key value":        "protocol=1\nunix=/tmp/a.sock\nauth\nEOF\n",
		"unknown key":          "protocol=1\nunix=/tmp/a.sock\nauth=" + testToken() + "\nextra=1\nEOF\n",
		"duplicate key":        "protocol=1\nprotocol=1\nunix=/tmp/a.sock\nauth=" + testToken() + "\nEOF\n",
		"unknown protocol":     "protocol=2\nunix=/tmp/a.sock\nauth=" + testToken() + "\nEOF\n",
		"non numeric protocol": "protocol=abc\nunix=/tmp/a.sock\nauth=" + testToken() + "\nEOF\n",
		"missing protocol":     "unix=/tmp/a.sock\nauth=" + testToken() + "\nEOF\n",
		"missing auth":         "protocol=1\nunix=/tmp/a.sock\nEOF\n",
		"missing endpoint":     "protocol=1\nauth=" + testToken() + "\nEOF\n",
		"both endpoints":       "protocol=1\nunix=/tmp/a.sock\nsock=12345\nauth=" + testToken() + "\nEOF\n",
		"relative unix path":   "protocol=1\nunix=tmp/a.sock\nauth=" + testToken() + "\nEOF\n",
		"oversize unix path":   "protocol=1\nunix=/" + strings.Repeat("a", maxUnixPathLen) + "\nauth=" + testToken() + "\nEOF\n",
		"port zero":            "protocol=1\nsock=0\nauth=" + testToken() + "\nEOF\n",
		"port range":           "protocol=1\nsock=65536\nauth=" + testToken() + "\nEOF\n",
		"port leading zero":    "protocol=1\nsock=01234\nauth=" + testToken() + "\nEOF\n",
		"port not digits":      "protocol=1\nsock=12a45\nauth=" + testToken() + "\nEOF\n",
		"token too short":      "protocol=1\nunix=/tmp/a.sock\nauth=QiQi\nEOF\n",
		"token std alphabet":   "protocol=1\nunix=/tmp/a.sock\nauth=" + strings.Repeat("+", 43) + "\nEOF\n",
		"token padded":         "protocol=1\nunix=/tmp/a.sock\nauth=" + strings.TrimSuffix(testToken(), "i") + "i=\nEOF\n",
	}
	for name, content := range cases {
		if _, err := Parse([]byte(content)); err == nil {
			t.Fatalf("%s: Parse unexpectedly succeeded", name)
		}
	}
}

func TestParseRejectsOversizeContent(t *testing.T) {
	content := "protocol=1\nunix=/tmp/" + strings.Repeat("a", MaxFileSize) + "\nauth=" + testToken() + "\nEOF\n"
	if _, err := Parse([]byte(content)); err == nil {
		t.Fatal("Parse unexpectedly succeeded for oversize content")
	}
}
