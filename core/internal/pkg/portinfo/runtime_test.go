package portinfo

import (
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"testing"
)

func TestDefaultRuntimeDirFormat(t *testing.T) {
	dir, err := DefaultRuntimeDir()
	if err != nil {
		t.Fatalf("DefaultRuntimeDir: %v", err)
	}
	base := filepath.Base(dir)
	if !strings.HasPrefix(base, RuntimeDirPrefix) {
		t.Fatalf("dir base %q missing prefix %q", base, RuntimeDirPrefix)
	}
	if runtime.GOOS != "windows" {
		want := filepath.Join("/tmp", RuntimeDirPrefix+strconv.Itoa(os.Getuid()))
		if dir != want {
			t.Fatalf("DefaultRuntimeDir = %q, want %q", dir, want)
		}
	}
}

func TestEnsureRuntimeDirCreatesAndReuses(t *testing.T) {
	root := t.TempDir()
	dir, err := runtimeDirInRoot(root)
	if err != nil {
		t.Fatalf("runtimeDirInRoot: %v", err)
	}

	if ensureErr := EnsureRuntimeDir(dir); ensureErr != nil {
		t.Fatalf("first EnsureRuntimeDir: %v", ensureErr)
	}
	if ensureErr := EnsureRuntimeDir(dir); ensureErr != nil {
		t.Fatalf("second EnsureRuntimeDir: %v", ensureErr)
	}

	if runtime.GOOS != "windows" {
		info, statErr := os.Stat(dir)
		if statErr != nil {
			t.Fatalf("Stat: %v", statErr)
		}
		if perm := info.Mode().Perm(); perm != RuntimeDirPerm {
			t.Fatalf("dir perm = %#o, want %#o", perm, RuntimeDirPerm)
		}
	}
}

func TestValidateRuntimeDirRejectsFileAndSymlink(t *testing.T) {
	root := t.TempDir()

	filePath := filepath.Join(root, "not-a-dir")
	if err := os.WriteFile(filePath, []byte("x"), 0o600); err != nil {
		t.Fatal(err)
	}
	if err := ValidateRuntimeDir(filePath); err == nil {
		t.Fatal("expected regular file to be rejected")
	}

	if runtime.GOOS != "windows" {
		realDir := filepath.Join(root, "real")
		if err := os.Mkdir(realDir, RuntimeDirPerm); err != nil {
			t.Fatal(err)
		}
		linkPath := filepath.Join(root, "link")
		if err := os.Symlink(realDir, linkPath); err != nil {
			t.Fatal(err)
		}
		if err := ValidateRuntimeDir(linkPath); err == nil {
			t.Fatal("expected symlink runtime dir to be rejected")
		}

		badPermDir := filepath.Join(root, "bad-perm")
		if err := os.Mkdir(badPermDir, 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.Chmod(badPermDir, 0o755); err != nil {
			t.Fatal(err)
		}
		if err := ValidateRuntimeDir(badPermDir); err == nil {
			t.Fatal("expected 0755 runtime dir to be rejected")
		}
	}
}

func TestValidateSocketPathLimit(t *testing.T) {
	if runtime.GOOS == "windows" {
		if err := ValidateSocketPath(strings.Repeat("a", 200)); err != nil {
			t.Fatalf("windows ValidateSocketPath should be no-op, got %v", err)
		}
		return
	}

	shortDir := "/tmp/slc-ok"
	if err := ValidateSocketPath(shortDir); err != nil {
		t.Fatalf("ValidateSocketPath(%q): %v", shortDir, err)
	}

	// core.sock 长度为 9，加上 "/" 共 10 字节，父目录达到 94 字节时总长 104 > 103。
	overDir := "/" + strings.Repeat("a", MaxUnixPathLen-len(CoreSocketFileName))
	if err := ValidateSocketPath(overDir); err == nil {
		t.Fatalf("expected socket path > %d bytes to be rejected", MaxUnixPathLen)
	}
}
