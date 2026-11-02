package watcher

import (
	"os"
	"path/filepath"
	"runtime"
	"testing"

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

	created, err := CreateSaveLinks([]*savev1.SaveRecord{save}, filesDir)
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

	created, err = CreateSaveLinks([]*savev1.SaveRecord{save}, filesDir)
	if err != nil {
		t.Fatal(err)
	}
	if created != 0 {
		t.Fatalf("second call created = %d, want 0", created)
	}
}

func TestCreateSaveLinksRequiresFilesDir(t *testing.T) {
	if _, err := CreateSaveLinks(nil, ""); err == nil {
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
