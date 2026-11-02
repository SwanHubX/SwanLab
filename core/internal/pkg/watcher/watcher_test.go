package watcher

import (
	"errors"
	"os"
	"path/filepath"
	"sync"
	"testing"
	"time"

	"go.uber.org/goleak"

	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

// fakeTimer 记录防抖回调，由测试显式触发。
type fakeTimer struct {
	fn   func()
	done bool
}

func (t *fakeTimer) Stop() bool {
	if t.done {
		return false
	}
	t.done = true
	return true
}

type fakeClock struct {
	mu     sync.Mutex
	delay  time.Duration
	timers []*fakeTimer
}

func (c *fakeClock) after(d time.Duration, fn func()) debounceTimer {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.delay = d
	timer := &fakeTimer{fn: fn}
	c.timers = append(c.timers, timer)
	return timer
}

func (c *fakeClock) latest() *fakeTimer {
	c.mu.Lock()
	defer c.mu.Unlock()
	return c.timers[len(c.timers)-1]
}

func (c *fakeClock) fire(timer *fakeTimer) {
	if timer.Stop() {
		timer.fn()
	}
}

func withTimerFactory(factory timerFactory) Option {
	return func(w *FileWatcher) { w.afterFunc = factory }
}

func newSave(name, source, target string, policy savev1.SavePolicy) *savev1.SaveRecord {
	return &savev1.SaveRecord{
		Name:       name,
		SourcePath: source,
		TargetPath: target,
		Policy:     policy,
	}
}

// addEntry 直接注册已存在的文件，绕过 fsnotify，便于确定性测试防抖逻辑。
func addEntry(t *testing.T, w *FileWatcher, name, source string, policy savev1.SavePolicy) string {
	t.Helper()
	info, ok := statFile(source)
	if !ok {
		t.Fatalf("stat file %s", source)
	}
	target := filepath.Join(w.filesDir, name)
	w.mu.Lock()
	defer w.mu.Unlock()
	w.entries[target] = &fileEntry{
		name:       name,
		sourcePath: source,
		targetPath: target,
		policy:     policy,
		info:       info,
	}
	return target
}

func TestStatFileMissing(t *testing.T) {
	if _, ok := statFile(filepath.Join(t.TempDir(), "missing.pt")); ok {
		t.Fatal("missing file must have no state")
	}
}

func TestSameFileTracksContent(t *testing.T) {
	path := filepath.Join(t.TempDir(), "model.pt")
	if err := os.WriteFile(path, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	first, ok := statFile(path)
	if !ok {
		t.Fatal("existing file must have state")
	}
	if err := os.WriteFile(path, []byte("v2-longer"), 0o644); err != nil {
		t.Fatal(err)
	}
	second, ok := statFile(path)
	if !ok {
		t.Fatal("existing file must have state")
	}
	if sameFile(first, second) {
		t.Fatal("content change was not detected")
	}
}

func TestSameFileDetectsAtomicReplace(t *testing.T) {
	dir := t.TempDir()
	source := filepath.Join(dir, "model.pt")
	if err := os.WriteFile(source, []byte("aa"), 0o644); err != nil {
		t.Fatal(err)
	}
	first, ok := statFile(source)
	if !ok {
		t.Fatal("existing file must have state")
	}

	temp := filepath.Join(dir, "model.pt.tmp")
	if err := os.WriteFile(temp, []byte("bb"), 0o644); err != nil {
		t.Fatal(err)
	}
	// 保持 mtime 与尺寸不变，仅文件身份变化。
	if err := os.Chtimes(temp, first.ModTime(), first.ModTime()); err != nil {
		t.Fatal(err)
	}
	if err := os.Rename(temp, source); err != nil {
		t.Fatal(err)
	}
	if err := os.Chtimes(source, first.ModTime(), first.ModTime()); err != nil {
		t.Fatal(err)
	}
	second, ok := statFile(source)
	if !ok {
		t.Fatal("existing file must have state")
	}
	if sameFile(first, second) {
		t.Fatal("atomic replacement was not detected")
	}
}

func TestNewValidatesArguments(t *testing.T) {
	callback := func(*savev1.SaveRecord) {}
	if _, err := New("", callback); err == nil {
		t.Fatal("empty files dir must fail")
	}
	if _, err := New(t.TempDir(), nil); err == nil {
		t.Fatal("nil callback must fail")
	}
	if _, err := New(t.TempDir(), callback, WithDebounce(0)); err == nil {
		t.Fatal("non-positive debounce must fail")
	}
}

func TestDebounceCoalescesChanges(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}

	clock := &fakeClock{}
	var mu sync.Mutex
	var records []*savev1.SaveRecord
	w, err := New(filepath.Join(root, "files"), func(record *savev1.SaveRecord) {
		mu.Lock()
		records = append(records, record)
		mu.Unlock()
	}, WithDebounce(2*time.Second), withTimerFactory(clock.after))
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()

	target := addEntry(t, w, "model.pt", source, savev1.SavePolicy_SAVE_POLICY_LIVE)
	w.notify(source)
	first := clock.latest()
	if clock.delay != 2*time.Second {
		t.Fatalf("debounce delay = %v, want 2s", clock.delay)
	}

	w.notify(source)
	second := clock.latest()
	if first == second || !first.done {
		t.Fatal("second event must reset the pending timer")
	}

	if err := os.WriteFile(source, []byte("v2"), 0o644); err != nil {
		t.Fatal(err)
	}
	clock.fire(second)
	mu.Lock()
	if len(records) != 1 {
		t.Fatalf("callbacks = %d, want 1", len(records))
	}
	record := records[0]
	mu.Unlock()
	if record.GetName() != "model.pt" || record.GetSourcePath() != source || record.GetTargetPath() != target {
		t.Fatalf("unexpected record: %v", record)
	}
	if record.GetPolicy() != savev1.SavePolicy_SAVE_POLICY_LIVE {
		t.Fatalf("policy = %v, want LIVE", record.GetPolicy())
	}

	// 文件状态未变化时重复触发不发回调。
	w.notify(source)
	clock.fire(clock.latest())
	mu.Lock()
	if len(records) != 1 {
		t.Fatalf("unchanged file triggered callback, total %d", len(records))
	}
	mu.Unlock()

	if err := os.WriteFile(source, []byte("v2-longer"), 0o644); err != nil {
		t.Fatal(err)
	}
	w.notify(source)
	clock.fire(clock.latest())
	mu.Lock()
	if len(records) != 2 {
		t.Fatalf("changed file callbacks = %d, want 2", len(records))
	}
	mu.Unlock()
}

func TestStaleTimerDoesNotClobberNewSchedule(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	calls := 0
	clock := &fakeClock{}
	w, err := New(filepath.Join(root, "files"), func(*savev1.SaveRecord) { calls++ }, withTimerFactory(clock.after))
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()

	addEntry(t, w, "model.pt", source, savev1.SavePolicy_SAVE_POLICY_LIVE)
	w.notify(source)
	stale := clock.latest()
	w.notify(source)
	latest := clock.latest()

	// 旧定时器已过期待执行，不能消耗新调度。
	stale.fn()
	if calls != 0 {
		t.Fatalf("stale timer triggered %d callbacks", calls)
	}
	if err := os.WriteFile(source, []byte("v2"), 0o644); err != nil {
		t.Fatal(err)
	}
	clock.fire(latest)
	if calls != 1 {
		t.Fatalf("new timer callbacks = %d, want 1", calls)
	}
}

func TestMissingFileKeepsRegistration(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	calls := 0
	clock := &fakeClock{}
	w, err := New(filepath.Join(root, "files"), func(*savev1.SaveRecord) { calls++ }, withTimerFactory(clock.after))
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()

	target := addEntry(t, w, "model.pt", source, savev1.SavePolicy_SAVE_POLICY_LIVE)
	if err := os.Remove(source); err != nil {
		t.Fatal(err)
	}
	w.notify(source)
	clock.fire(clock.latest())
	if calls != 0 {
		t.Fatalf("missing file triggered %d callbacks", calls)
	}
	w.mu.Lock()
	entry, ok := w.entries[target]
	if !ok || entry.info != nil {
		w.mu.Unlock()
		t.Fatalf("missing file must keep registration and clear state: %+v", entry)
	}
	w.mu.Unlock()

	if err := os.WriteFile(source, []byte("recreated"), 0o644); err != nil {
		t.Fatal(err)
	}
	w.notify(source)
	clock.fire(clock.latest())
	if calls != 1 {
		t.Fatalf("recreated file callbacks = %d, want 1", calls)
	}
}

func TestStopFlushesPendingChanges(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	calls := 0
	clock := &fakeClock{}
	w, err := New(filepath.Join(root, "files"), func(*savev1.SaveRecord) { calls++ },
		WithDebounce(time.Hour), withTimerFactory(clock.after))
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()

	addEntry(t, w, "model.pt", source, savev1.SavePolicy_SAVE_POLICY_LIVE)
	w.notify(source)
	pending := clock.latest()
	if err := os.WriteFile(source, []byte("v2"), 0o644); err != nil {
		t.Fatal(err)
	}

	w.Stop()
	if calls != 1 {
		t.Fatalf("stop callbacks = %d, want 1", calls)
	}
	if !pending.done {
		t.Fatal("stop must cancel pending timer")
	}

	w.notify(source)
	if clock.latest() != pending {
		t.Fatal("stopped watcher must not schedule new timers")
	}
	w.Stop()
	if calls != 1 {
		t.Fatalf("second stop callbacks = %d, want 1", calls)
	}
}

func TestStopChecksAllRegisteredFiles(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	calls := 0
	clock := &fakeClock{}
	w, err := New(filepath.Join(root, "files"), func(*savev1.SaveRecord) { calls++ },
		WithDebounce(time.Hour), withTimerFactory(clock.after))
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()

	addEntry(t, w, "model.pt", source, savev1.SavePolicy_SAVE_POLICY_LIVE)
	if err := os.WriteFile(source, []byte("v2"), 0o644); err != nil {
		t.Fatal(err)
	}

	// 没有 notify/timer 也要在 Stop 时检查到最终状态。
	w.Stop()
	if calls != 1 {
		t.Fatalf("stop callbacks = %d, want 1", calls)
	}
}

func TestRegisterLiveSkipsOtherPolicies(t *testing.T) {
	w, err := New(t.TempDir(), func(*savev1.SaveRecord) {})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()

	saves := []*savev1.SaveRecord{
		newSave("now.pt", "/tmp/now.pt", "", savev1.SavePolicy_SAVE_POLICY_NOW),
		newSave("end.pt", "/tmp/end.pt", "", savev1.SavePolicy_SAVE_POLICY_END),
	}
	if err := w.RegisterLive(saves); err != nil {
		t.Fatal(err)
	}
	w.mu.Lock()
	defer w.mu.Unlock()
	if w.fsWatcher != nil || len(w.entries) != 0 {
		t.Fatal("non-LIVE records must not start the watcher")
	}
}

func TestRegisterLiveRequiresExistingSourceDir(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "missing", "model.pt")
	w, err := New(filepath.Join(root, "files"), func(*savev1.SaveRecord) {})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Stop()

	save := newSave("model.pt", source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)
	if err := w.RegisterLive([]*savev1.SaveRecord{save}); err == nil {
		t.Fatal("watching a missing parent dir must fail")
	}
	w.mu.Lock()
	defer w.mu.Unlock()
	if len(w.entries) != 0 {
		t.Fatalf("failed registration left %d entries", len(w.entries))
	}
}

func TestRegisterLiveAfterStop(t *testing.T) {
	root := t.TempDir()
	source := filepath.Join(root, "model.pt")
	if err := os.WriteFile(source, []byte("v1"), 0o644); err != nil {
		t.Fatal(err)
	}
	w, err := New(filepath.Join(root, "files"), func(*savev1.SaveRecord) {})
	if err != nil {
		t.Fatal(err)
	}
	w.Stop()

	save := newSave("model.pt", source, "", savev1.SavePolicy_SAVE_POLICY_LIVE)
	if err := w.RegisterLive([]*savev1.SaveRecord{save}); !errors.Is(err, errStopped) {
		t.Fatalf("error = %v, want errStopped", err)
	}
}

func TestNullWatcherIsNoop(t *testing.T) {
	w := NullWatcher{}
	save := newSave("model.pt", "/tmp/model.pt", "", savev1.SavePolicy_SAVE_POLICY_LIVE)
	if err := w.RegisterLive([]*savev1.SaveRecord{save}); err != nil {
		t.Fatal(err)
	}
	w.Stop()
}
