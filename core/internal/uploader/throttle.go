package uploader

import (
	"sync"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

const (
	// DefaultWarningInterval 是失败告警的最小间隔。
	DefaultWarningInterval = 30 * time.Second
	// uploadFailureWarning 是失败告警文案。
	uploadFailureWarning = "Upload failed due to network or server issues, retrying automatically."
	// uploadRecoveredNotice 是恢复提示文案。
	uploadRecoveredNotice = "Upload recovered, resuming normally."
)

// WarningThrottle 控制失败告警频率：首次立即，其后每个间隔至多一条；恢复时提示一次。
type WarningThrottle struct {
	mu       sync.Mutex
	now      func() time.Time
	interval time.Duration
	warnf    func(string)
	infof    func(string)

	lastWarnAt time.Time
	inFailure  bool
}

// NewWarningThrottle 构造默认节流器（间隔 30s）。
func NewWarningThrottle() *WarningThrottle {
	return newWarningThrottle(time.Now, DefaultWarningInterval)
}

// newWarningThrottle 允许注入时钟与间隔。
func newWarningThrottle(now func() time.Time, interval time.Duration) *WarningThrottle {
	return &WarningThrottle{
		now:        now,
		interval:   interval,
		warnf:      func(msg string) { console.Warning(msg) },
		infof:      func(msg string) { console.Info(msg) },
		lastWarnAt: now().Add(-interval), // 使首次调用必定触发
	}
}

// Warn 记录一次失败，必要时输出告警。
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

// Reset 标记已恢复：故障期内提示一次，并复位告警窗口。
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
