package watcher

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/pkg/fs"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func TestCreateSaveLinks(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("weights"), 0o644); err != nil {
		t.Fatal(err)
	}
	filesDir := filepath.Join(root, "run", "files")
	save := newSave("nested/model.pt", source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)

	created, err := CreateSaveLinks(context.Background(), []*savev1.SaveRecord{save}, filesDir, 0)
	if err != nil {
		t.Fatal(err)
	}
	if created != 1 {
		t.Fatalf("created = %d, want 1", created)
	}
	target := filepath.Join(filesDir, "nested", "model.pt")
	if save.GetTargetPath() != target {
		t.Fatalf("target path = %q, want %q", save.GetTargetPath(), target)
	}
	info, err := os.Lstat(target)
	if err != nil {
		t.Fatal(err)
	}
	if runtime.GOOS != "windows" && info.Mode()&os.ModeSymlink == 0 {
		t.Fatalf("target mode = %v, want symlink", info.Mode())
	}
	content, err := os.ReadFile(target)
	if err != nil {
		t.Fatal(err)
	}
	if string(content) != "weights" {
		t.Fatalf("target content = %q", content)
	}

	created, err = CreateSaveLinks(context.Background(), []*savev1.SaveRecord{save}, filesDir, 0)
	if err != nil {
		t.Fatal(err)
	}
	if created != 0 {
		t.Fatalf("second call created = %d, want 0", created)
	}
}

func TestCreateSaveLinksUsesWindowsParts(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("weights"), 0o644); err != nil {
		t.Fatal(err)
	}
	filesDir := filepath.Join(root, "run", "files")
	save := newSave(`nested\model.pt`, source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)

	if _, err := CreateSaveLinks(context.Background(), []*savev1.SaveRecord{save}, filesDir, 0); err != nil {
		t.Fatal(err)
	}
	want := filepath.Join(filesDir, "nested", "model.pt")
	if save.GetTargetPath() != want {
		t.Fatalf("target path = %q, want %q", save.GetTargetPath(), want)
	}
	if got := fs.MirrorPath(filesDir, save.GetName()); got != want {
		t.Fatalf("sender fallback path = %q, want %q", got, want)
	}
	if _, err := os.Stat(want); err != nil {
		t.Fatalf("mirror missing: %v", err)
	}
}

func TestCreateSaveLinksRequiresFilesDir(t *testing.T) {
	if _, err := CreateSaveLinks(context.Background(), nil, "", 0); err == nil {
		t.Fatal("empty files dir must fail")
	}
}

func TestResolveSourcePrefersSourceForRegularTarget(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "source.pt")
	target := filepath.Join(root, "copy.pt")
	for _, path := range []string{source, target} {
		if err := os.WriteFile(path, []byte("v1"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if got := resolveSource(target, source); got != source {
		t.Fatalf("resolved source = %q, want %q", got, source)
	}
	expected, err := filepath.EvalSymlinks(target)
	if err != nil {
		t.Fatal(err)
	}
	if got := resolveSource(target, ""); got != expected {
		t.Fatalf("resolved source without fallback = %q, want %q", got, expected)
	}
}

func TestResolveSourceFollowsSymlink(t *testing.T) {
	root := t.TempDir()
	realSource := filepath.Join(root, "real.pt")
	if err := os.WriteFile(realSource, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	target := filepath.Join(root, "link.pt")
	if err := os.Symlink(realSource, target); err != nil {
		t.Skipf("symlink not supported: %v", err)
	}
	expected, err := filepath.EvalSymlinks(realSource)
	if err != nil {
		t.Fatal(err)
	}
	if got := resolveSource(target, ""); got != expected {
		t.Fatalf("resolved source = %q, want %q", got, expected)
	}
}

// 回退拷贝前就要按保存上限拒绝，不能白拷贝一遍再被上传端丢弃。
func TestCopyFileRejectsOversizedFile(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "big.bin")
	if err := os.WriteFile(source, make([]byte, 2048), 0o644); err != nil {
		t.Fatal(err)
	}
	target := filepath.Join(root, "copy.bin")

	err := copyFile(context.Background(), source, target, 1024)

	if err == nil || !strings.Contains(err.Error(), "exceeds the save size limit") {
		t.Fatalf("copyFile err = %v, want size limit error", err)
	}
	if _, statErr := os.Stat(target); !os.IsNotExist(statErr) {
		t.Fatalf("oversized file must not be copied, stat err = %v", statErr)
	}
}

// 拷贝过程中取消要立刻中止，且不能留下半截文件。
func TestCopyFileStopsOnCanceledContext(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.bin")
	if err := os.WriteFile(source, make([]byte, 4<<20), 0o644); err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	target := filepath.Join(root, "copy.bin")

	if err := copyFile(ctx, source, target, 0); !errors.Is(err, context.Canceled) {
		t.Fatalf("copyFile err = %v, want context.Canceled", err)
	}
	if _, statErr := os.Stat(target); !os.IsNotExist(statErr) {
		t.Fatalf("canceled copy must not leave a partial file, stat err = %v", statErr)
	}
}

func TestCreateSaveLinksStopsOnCanceledContext(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("weights"), 0o644); err != nil {
		t.Fatal(err)
	}
	save := newSave("model.pt", source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)
	ctx, cancel := context.WithCancel(context.Background())
	cancel()

	if _, err := CreateSaveLinks(ctx, []*savev1.SaveRecord{save}, filepath.Join(root, "files"), 0); !errors.Is(
		err, context.Canceled,
	) {
		t.Fatalf("CreateSaveLinks err = %v, want context.Canceled", err)
	}
}
