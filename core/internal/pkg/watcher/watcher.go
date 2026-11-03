// Package watcher 监听 LIVE 策略的保存文件，在文件稳定后通知调用方重新上传。
//
// 监听软链接指向的真实源文件父目录，以捕获外部编辑与原子替换；按 mtime、
// 大小与文件身份去重，每个文件独立防抖，Stop 会同步刷出待发变更。
//
// 本包只依赖标准库、fsnotify、save proto 与 internal/pkg/fs，不依赖业务包。
package watcher

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"sync"
	"time"

	"github.com/fsnotify/fsnotify"

	"github.com/swanhubx/swanlab/core/internal/pkg/fs"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

const (
	// DefaultDebounce 是文件停止变化后等待再次上传的时间。
	DefaultDebounce = 1500 * time.Millisecond
	// fsEventMask 是触发重新上传的事件集合。
	fsEventMask = fsnotify.Create | fsnotify.Write | fsnotify.Chmod | fsnotify.Rename
)

var errStopped = errors.New("watcher: stopped")

// Watcher 是 LIVE 文件监听的对外接口，FileWatcher 与 NullWatcher 均实现它。
type Watcher interface {
	// RegisterLive 为 policy=SAVE_POLICY_LIVE 的记录注册监听，重复注册保持不变。
	RegisterLive(saves []*savev1.SaveRecord) error
	// Stop 停止监听并同步刷出尚在防抖窗口内的变更。
	Stop()
}

// FileWatcher 基于 fsnotify 监听真实源文件所在目录。
// 原子替换会更换 inode，直接监听文件会丢失后续事件。
type FileWatcher struct {
	filesDir  string
	onChange  func(*savev1.SaveRecord)
	debounce  time.Duration
	afterFunc timerFactory

	mu          sync.Mutex
	cond        *sync.Cond
	stopped     bool
	entries     map[string]*fileEntry // target path → entry
	timers      map[string]scheduled
	nextGen     uint64
	watchedDirs map[string]struct{}
	active      int
	fsWatcher   *fsnotify.Watcher
	loopWG      sync.WaitGroup
}

// fileEntry 是单个 LIVE 文件的注册信息。
type fileEntry struct {
	name       string
	sourcePath string
	targetPath string
	policy     savev1.SavePolicy
	info       os.FileInfo // nil 表示注册时文件缺失
}

// statFile 通过文件句柄取状态；Windows 下文件索引在取状态时即固定。
// 文件缺失时 ok 为 false。
func statFile(path string) (os.FileInfo, bool) {
	file, err := os.Open(path)
	if err != nil {
		return nil, false
	}
	defer func() { _ = file.Close() }()

	info, err := file.Stat()
	if err != nil {
		return nil, false
	}
	return info, true
}

// sameFile 判断文件身份、大小与 mtime 是否都未变化。
func sameFile(previous, current os.FileInfo) bool {
	return previous.Size() == current.Size() &&
		previous.ModTime().Equal(current.ModTime()) &&
		os.SameFile(previous, current)
}

// Option 调整 FileWatcher 的行为。
type Option func(*FileWatcher)

// WithDebounce 设置尾部防抖时长，非正值会被 New 拒绝。
func WithDebounce(d time.Duration) Option {
	return func(w *FileWatcher) { w.debounce = d }
}

type debounceTimer interface {
	Stop() bool
}

// scheduled 用自增代号标记一次调度，防止已过期但仍在执行的定时器影响新调度。
type scheduled struct {
	timer debounceTimer
	gen   uint64
}

type timerFactory func(time.Duration, func()) debounceTimer

func newTimer(d time.Duration, fn func()) debounceTimer { return time.AfterFunc(d, fn) }

// New 创建文件监听器。filesDir 是 run 目录下的 files 镜像目录。
// 监听器惰性启动：首次收到 LIVE 记录时才创建 fsnotify。
func New(filesDir string, onChange func(*savev1.SaveRecord), opts ...Option) (*FileWatcher, error) {
	if filesDir == "" {
		return nil, errors.New("watcher: files dir is required")
	}
	if onChange == nil {
		return nil, errors.New("watcher: on-change callback is required")
	}
	abs, err := filepath.Abs(filesDir)
	if err != nil {
		return nil, fmt.Errorf("watcher: resolve files dir: %w", err)
	}
	w := &FileWatcher{
		filesDir:    abs,
		onChange:    onChange,
		debounce:    DefaultDebounce,
		afterFunc:   newTimer,
		entries:     make(map[string]*fileEntry),
		timers:      make(map[string]scheduled),
		watchedDirs: make(map[string]struct{}),
	}
	w.cond = sync.NewCond(&w.mu)
	for _, opt := range opts {
		opt(w)
	}
	if w.debounce <= 0 {
		return nil, fmt.Errorf("watcher: debounce must be positive, got %v", w.debounce)
	}
	return w, nil
}

// RegisterLive 为 LIVE 记录注册监听。非 LIVE 记录会被忽略。
// 批量注册按事务处理：任一条失败时撤销本批新增的条目与目录，既有监听保持不变。
func (w *FileWatcher) RegisterLive(saves []*savev1.SaveRecord) error {
	w.mu.Lock()
	defer w.mu.Unlock()

	var (
		addedTargets []string
		addedDirs    []string
		errs         []error
	)
	for _, save := range saves {
		if save.GetPolicy() != savev1.SavePolicy_SAVE_POLICY_LIVE {
			continue
		}
		target, dir, added, err := w.registerLocked(save)
		if added {
			addedTargets = append(addedTargets, target)
		}
		if dir != "" {
			addedDirs = append(addedDirs, dir)
		}
		if err != nil {
			errs = append(errs, fmt.Errorf("watch %q: %w", save.GetName(), err))
			break
		}
	}
	if len(errs) > 0 {
		w.rollbackLocked(addedTargets, addedDirs)
		return errors.Join(errs...)
	}
	return nil
}

// registerLocked 注册单个文件，返回新增的条目与目录供失败回滚；调用方需持有 mu。
func (w *FileWatcher) registerLocked(save *savev1.SaveRecord) (target, newDir string, added bool, err error) {
	target = save.GetTargetPath()
	if target == "" {
		target = fs.MirrorPath(w.filesDir, save.GetName())
	}
	target = filepath.Clean(target)
	source := resolveSource(target, save.GetSourcePath())

	if w.stopped {
		return target, "", false, errStopped
	}
	if _, ok := w.entries[target]; ok {
		return target, "", false, nil
	}
	if err := w.startLocked(); err != nil {
		return target, "", false, err
	}
	dir := filepath.Dir(source)
	if _, ok := w.watchedDirs[dir]; !ok {
		if err := w.fsWatcher.Add(dir); err != nil {
			return target, "", false, fmt.Errorf("watch dir %q: %w", dir, err)
		}
		w.watchedDirs[dir] = struct{}{}
		newDir = dir
	}
	info, _ := statFile(source)
	w.entries[target] = &fileEntry{
		name:       save.GetName(),
		sourcePath: source,
		targetPath: target,
		policy:     save.GetPolicy(),
		info:       info,
	}
	return target, newDir, true, nil
}

// rollbackLocked 撤销本批新增的监听，不影响此前已注册的文件；调用方需持有 mu。
func (w *FileWatcher) rollbackLocked(targets, dirs []string) {
	for _, target := range targets {
		delete(w.entries, target)
		if timer, ok := w.timers[target]; ok {
			timer.timer.Stop()
			delete(w.timers, target)
		}
	}
	for _, dir := range dirs {
		if w.fsWatcher != nil {
			_ = w.fsWatcher.Remove(dir)
		}
		delete(w.watchedDirs, dir)
	}
}

// resolveSource 解析软链接指向的真实文件；目标不是软链接时优先用源路径，
// 这样软链接创建失败、退化为复制的场景仍能监听原始文件。
func resolveSource(target, fallback string) string {
	if info, err := os.Lstat(target); err == nil && info.Mode()&os.ModeSymlink == 0 && fallback != "" {
		return filepath.Clean(fallback)
	}
	if resolved, err := filepath.EvalSymlinks(target); err == nil {
		return resolved
	}
	if fallback != "" {
		return filepath.Clean(fallback)
	}
	return target
}

// startLocked 在第一次注册时创建 fsnotify 并启动事件循环。
func (w *FileWatcher) startLocked() error {
	if w.fsWatcher != nil {
		return nil
	}
	fsWatcher, err := fsnotify.NewWatcher()
	if err != nil {
		return fmt.Errorf("create fsnotify watcher: %w", err)
	}
	w.fsWatcher = fsWatcher
	w.loopWG.Add(1)
	go w.loop(fsWatcher)
	return nil
}

func (w *FileWatcher) loop(fsWatcher *fsnotify.Watcher) {
	defer w.loopWG.Done()
	for {
		select {
		case event, ok := <-fsWatcher.Events:
			if !ok {
				return
			}
			if event.Op&fsEventMask == 0 {
				continue
			}
			w.notify(filepath.Clean(event.Name))
		case _, ok := <-fsWatcher.Errors:
			// 底层监听错误不终止循环，避免单个目录的异常影响其它文件。
			if !ok {
				return
			}
		}
	}
}

// notify 为源路径命中的文件重置防抖定时器。
func (w *FileWatcher) notify(source string) {
	w.mu.Lock()
	defer w.mu.Unlock()
	if w.stopped {
		return
	}
	for target, entry := range w.entries {
		if entry.sourcePath == source {
			w.scheduleLocked(target)
		}
	}
}

func (w *FileWatcher) scheduleLocked(target string) {
	if old, ok := w.timers[target]; ok {
		old.timer.Stop()
	}
	w.nextGen++
	gen := w.nextGen
	timer := w.afterFunc(w.debounce, func() { w.fire(target, gen) })
	w.timers[target] = scheduled{timer: timer, gen: gen}
}

// fire 在防抖到期后比较文件状态，变化时触发回调。
func (w *FileWatcher) fire(target string, gen uint64) {
	w.mu.Lock()
	if w.stopped {
		w.mu.Unlock()
		return
	}
	current, ok := w.timers[target]
	if !ok || current.gen != gen {
		w.mu.Unlock()
		return
	}
	delete(w.timers, target)
	record := w.takeChangeLocked(target)
	if record == nil {
		w.mu.Unlock()
		return
	}
	w.active++
	w.mu.Unlock()

	w.onChange(record)

	w.mu.Lock()
	w.active--
	w.cond.Broadcast()
	w.mu.Unlock()
}

// takeChangeLocked 在文件状态变化时更新状态并返回待通知记录；调用方需持有 mu。
func (w *FileWatcher) takeChangeLocked(target string) *savev1.SaveRecord {
	entry := w.entries[target]
	if entry == nil {
		return nil
	}
	info, ok := statFile(entry.sourcePath)
	if !ok {
		// 文件缺失时保留注册并清空状态，重建后仍能触发。
		entry.info = nil
		return nil
	}
	if entry.info != nil && sameFile(entry.info, info) {
		return nil
	}
	entry.info = info
	return entry.record()
}

// Stop 停止监听，并在返回前同步刷出所有待发变更。
//
// 调用方不得持有 onChange 回调所需的锁：Stop 会等待在途回调，并同步执行收尾回调。
// 同理，不要在 onChange 内调用 Stop。
func (w *FileWatcher) Stop() {
	w.mu.Lock()
	if w.stopped {
		w.mu.Unlock()
		return
	}
	w.stopped = true
	for _, timer := range w.timers {
		timer.timer.Stop()
	}
	clear(w.timers)
	// 停止后对所有注册文件做一次状态检查，覆盖尚未进入 notify 的事件。
	targets := make([]string, 0, len(w.entries))
	for target := range w.entries {
		targets = append(targets, target)
	}
	fsWatcher := w.fsWatcher
	w.mu.Unlock()

	if fsWatcher != nil {
		_ = fsWatcher.Close()
		w.loopWG.Wait()
	}

	// 等待已进入回调的变更完成，避免收尾记录晚于 intake 封闭。
	w.mu.Lock()
	for w.active > 0 {
		w.cond.Wait()
	}
	w.mu.Unlock()

	sort.Strings(targets)
	for _, target := range targets {
		w.flush(target)
	}
}

// flush 在停止后同步处理一个待发目标。
func (w *FileWatcher) flush(target string) {
	w.mu.Lock()
	record := w.takeChangeLocked(target)
	w.mu.Unlock()
	if record != nil {
		w.onChange(record)
	}
}

func (e *fileEntry) record() *savev1.SaveRecord {
	return &savev1.SaveRecord{
		Name:       e.name,
		SourcePath: e.sourcePath,
		TargetPath: e.targetPath,
		Policy:     e.policy,
	}
}

// NullWatcher 用于 skip_store：不创建目录、不启动监听、不产生协程。
type NullWatcher struct{}

// RegisterLive 是空实现。
func (NullWatcher) RegisterLive([]*savev1.SaveRecord) error { return nil }

// Stop 是空实现。
func (NullWatcher) Stop() {}

var (
	_ Watcher = (*FileWatcher)(nil)
	_ Watcher = NullWatcher{}
)
