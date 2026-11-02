package watcher

import (
	"os"
	"path/filepath"
	"testing"
	"time"

	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func TestRegisterLiveIsIdempotent(t *testing.T) {
	root := t.TempDir()
	sourceDir := filepath.Join(root, "src")
	if err := os.MkdirAll(sourceDir, 0o755); err != nil {
		t.Fatal(err)
	}
	filesDir := filepath.Join(root, "run", "files")
	saves := make([]*savev1.SaveRecord, 0, 2)
	for _, name := range []string{"a.pt", "b.pt"} {
		source := filepath.Join(sourceDir, name)
		if err := os.WriteFile(source, []byte(name), 0o644); err != nil {
			t.Fatal(err)
		}
		saves = append(saves, newSave(name, source, "", savev1.SavePolicy_SAVE_POLICY_LIVE))
	}
	if _, err := CreateSaveLinks(saves, filesDir); err != nil {
		t.Fatal(err)
	}

	w, err := New(filesDir, func(*savev1.SaveRecord) {})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()
	if err := w.RegisterLive(saves); err != nil {
		t.Fatal(err)
	}
	if err := w.RegisterLive(saves); err != nil {
		t.Fatal(err)
	}

	w.mu.Lock()
	defer w.mu.Unlock()
	if len(w.entries) != 2 {
		t.Fatalf("entries = %d, want 2", len(w.entries))
	}
	if len(w.watchedDirs) != 1 {
		t.Fatalf("watched dirs = %d, want 1", len(w.watchedDirs))
	}
}

func TestWatchFileWrite(t *testing.T) {
	root := t.TempDir()
	sourceDir := filepath.Join(root, "src")
	if err := os.MkdirAll(sourceDir, 0o755); err != nil {
		t.Fatal(err)
	}
	source := filepath.Join(sourceDir, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	filesDir := filepath.Join(root, "run", "files")
	save := newSave("model.pt", source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)
	if _, err := CreateSaveLinks([]*savev1.SaveRecord{save}, filesDir); err != nil {
		t.Fatal(err)
	}

	records := make(chan *savev1.SaveRecord, 4)
	w, err := New(filesDir, func(record *savev1.SaveRecord) { records <- record }, WithDebounce(50*time.Millisecond))
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(w.Stop)
	if err = w.RegisterLive([]*savev1.SaveRecord{save}); err != nil {
		t.Fatal(err)
	}

	if err = os.WriteFile(source, []byte("v2"), 0o644); err != nil {
		t.Fatal(err)
	}
	expectedSource, err := filepath.EvalSymlinks(source)
	if err != nil {
		t.Fatal(err)
	}
	select {
	case record := <-records:
		if record.GetName() != "model.pt" || record.GetSourcePath() != expectedSource {
			t.Fatalf("unexpected record: %v", record)
		}
		if record.GetTargetPath() != save.GetTargetPath() {
			t.Fatalf("target = %q, want %q", record.GetTargetPath(), save.GetTargetPath())
		}
	case <-time.After(5 * time.Second):
		t.Fatal("timed out waiting for file write callback")
	}
}

func TestWatchAtomicReplace(t *testing.T) {
	root := t.TempDir()
	sourceDir := filepath.Join(root, "src")
	if err := os.MkdirAll(sourceDir, 0o755); err != nil {
		t.Fatal(err)
	}
	source := filepath.Join(sourceDir, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	filesDir := filepath.Join(root, "run", "files")
	save := newSave("model.pt", source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)
	if _, err := CreateSaveLinks([]*savev1.SaveRecord{save}, filesDir); err != nil {
		t.Fatal(err)
	}

	records := make(chan *savev1.SaveRecord, 4)
	w, err := New(filesDir, func(record *savev1.SaveRecord) { records <- record }, WithDebounce(50*time.Millisecond))
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(w.Stop)
	if err = w.RegisterLive([]*savev1.SaveRecord{save}); err != nil {
		t.Fatal(err)
	}

	temp := filepath.Join(sourceDir, "model.pt.tmp")
	if err = os.WriteFile(temp, []byte("v2-atomic"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err = os.Rename(temp, source); err != nil {
		t.Fatal(err)
	}
	select {
	case record := <-records:
		if record.GetName() != "model.pt" {
			t.Fatalf("unexpected record: %v", record)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("timed out waiting for atomic replace callback")
	}
}

func TestStopFlushesRealEvent(t *testing.T) {
	root := t.TempDir()
	sourceDir := filepath.Join(root, "src")
	if err := os.MkdirAll(sourceDir, 0o755); err != nil {
		t.Fatal(err)
	}
	source := filepath.Join(sourceDir, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	filesDir := filepath.Join(root, "run", "files")
	save := newSave("model.pt", source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)
	if _, err := CreateSaveLinks([]*savev1.SaveRecord{save}, filesDir); err != nil {
		t.Fatal(err)
	}

	records := make(chan *savev1.SaveRecord, 4)
	w, err := New(filesDir, func(record *savev1.SaveRecord) { records <- record }, WithDebounce(time.Hour))
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()
	scheduled := make(chan string, 4)
	w.scheduled = scheduled
	if err = w.RegisterLive([]*savev1.SaveRecord{save}); err != nil {
		t.Fatal(err)
	}

	if err = os.WriteFile(source, []byte("v2"), 0o644); err != nil {
		t.Fatal(err)
	}
	select {
	case <-scheduled:
	case <-time.After(5 * time.Second):
		t.Fatal("timed out waiting for debounce scheduling")
	}

	w.Stop()
	select {
	case <-records:
	default:
		t.Fatal("stop must flush the pending change before returning")
	}
}
