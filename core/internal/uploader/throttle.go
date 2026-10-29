package uploader

import (
	"sync"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

const (
	// DefaultWarningInterval 是同一次故障期内两条失败告警的最小间隔。
	DefaultWarningInterval = 30 * time.Second
	// uploadFailureWarning 对齐 Python UploadWarningThrottle.warn 的默认文案。
	uploadFailureWarning = "Upload failed due to network or server issues, retrying automatically."
	// uploadRecoveredNotice 对齐 Python UploadWarningThrottle.reset 的默认文案。
	uploadRecoveredNotice = "Upload recovered, resuming normally."
)

// WarningThrottle 控制上传失败告警的输出频率，并在恢复时提示一次。
//
// 首次失败立即告警，其后同一次故障期内每 DefaultWarningInterval 至多一条；
// 恢复提示只在故障期内输出过警报后触发，避免正常路径刷屏。
type WarningThrottle struct {
	mu       sync.Mutex
	now      func() time.Time
	interval time.Duration
	warnf    func(string)
	infof    func(string)

	lastWarnAt time.Time
	inFailure  bool
}

// NewWarningThrottle 构造默认节流器：间隔 30s，时钟取 time.Now。
func NewWarningThrottle() *WarningThrottle {
	return newWarningThrottle(time.Now, DefaultWarningInterval)
}

// newWarningThrottle 是测试入口，允许注入时钟与间隔。
func newWarningThrottle(now func() time.Time, interval time.Duration) *WarningThrottle {
	return &WarningThrottle{
		now:        now,
		interval:   interval,
		warnf:      func(msg string) { console.Warning(msg) },
		infof:      func(msg string) { console.Info(msg) },
		lastWarnAt: now().Add(-interval), // 使首次调用必定触发
	}
}

// Warn 记录一次上传失败，必要时输出告警。
func (t *WarningThrottle) Warn() {
	t.mu.Lock()
	t.inFailure = true
	now := t.now()
	warn := !now.Before(t.lastWarnAt.Add(t.interval))
	if warn {
		t.lastWarnAt = now
	}
	t.mu.Unlock()

	if warn {
		t.warnf(uploadFailureWarning)
	}
}

// Reset 标记上传已恢复：故障期内提示一次恢复，并复位告警窗口。
func (t *WarningThrottle) Reset() {
	t.mu.Lock()
	inFailure := t.inFailure
	t.inFailure = false
	t.lastWarnAt = t.now().Add(-t.interval)
	t.mu.Unlock()

	if inFailure {
		t.infof(uploadRecoveredNotice)
	}
}
