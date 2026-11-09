package manager

import (
	"context"
	"errors"
	"fmt"
	"sync"

	"github.com/google/uuid"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/cron/health"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

var (
	// ErrEmptySessionID 表示未提供有效的 session_id。
	ErrEmptySessionID = errors.New("session_id must not be empty")
	// ErrSessionNotFound 表示请求的 session_id 不存在或已被确认释放。
	ErrSessionNotFound = errors.New("unknown or released session_id")
	// ErrRunNotFinished 表示会话尚未排空，不能提前确认释放。
	ErrRunNotFinished = errors.New("run has not been drained for this handle")
	// ErrRunFailed 表示会话的收尾序列失败，不能按成功路径释放。
	ErrRunFailed = errors.New("run finish failed")
	// ErrNotRunning 表示会话已停止接受公共提交（排空中、已排空或收尾失败）。
	ErrNotRunning = errors.New("session is not accepting records")
	// ErrOverloaded 表示入口通道或上传配额在等待上限内没有腾出空间，可重试。
	ErrOverloaded = errors.New("session intake is overloaded")
	// ErrInvalidRecord 表示记录结构不合法（nil、空 key、未知 save 类型等）。
	ErrInvalidRecord = errors.New("invalid record")
	// ErrRequestTooLarge 表示单次 RPC 的记录数超过 record_batch。
	// 日志与 CUSTOM Save 无法跨块去重，必须在接纳前整体拒绝，重试不会成功。
	ErrRequestTooLarge = errors.New("request exceeds record_batch limit")
	// ErrTooManySessions 表示会话数已达上限，需要先释放已结束的会话。
	ErrTooManySessions = errors.New("session registry is full")
	// ErrInvalidSettings 表示启动配置不合法，属于调用方参数错误。
	ErrInvalidSettings = errors.New("invalid run settings")
	// ErrInvalidStartRecord 表示启动记录本身不合法（nil 或未知 resume 模式）。
	// 云端返回的元数据不完整属于协议错误，由 client.ErrInvalidResponse 表示。
	ErrInvalidStartRecord = errors.New("invalid start record")
	// ErrInvalidFinishRecord 表示未提供结束记录，属于调用方参数错误。
	ErrInvalidFinishRecord = errors.New("finish record must not be nil")
	// ErrControlUnreachable 表示启动预检无法访问控制面，快速失败避免认证重试堆积。
	ErrControlUnreachable = errors.New("control plane is unreachable")
)

const (
	// maxSessions 是单个 core 同时存在的会话上限。取值远高于单机多 rank 的进程数，
	// 触顶通常意味着调用方没有释放会话（例如循环 init 未 finish），此时明确失败
	// 而不是让注册表无限增长。
	maxSessions = 512

	// releasedCapacity 是已释放 session_id 的 FIFO 缓存上限，用于幂等 Confirm 重发。
	releasedCapacity = 128
)

// releasedCache 记录已成功释放的 session_id，按插入顺序淘汰最旧条目。
type releasedCache struct {
	entries map[string]struct{}
	order   []string
	next    int
}

func newReleasedCache() releasedCache {
	return releasedCache{entries: make(map[string]struct{}, releasedCapacity)}
}

func (c *releasedCache) has(id string) bool {
	_, ok := c.entries[id]
	return ok
}

func (c *releasedCache) add(id string) {
	if c.has(id) {
		return
	}
	if len(c.order) < releasedCapacity {
		c.order = append(c.order, id)
		c.entries[id] = struct{}{}
		return
	}
	delete(c.entries, c.order[c.next])
	c.order[c.next] = id
	c.entries[id] = struct{}{}
	c.next = (c.next + 1) % releasedCapacity
}

// Options 是 Manager 的构造参数。
type Options struct {
	// HTTP 配置共享连接池、重试与鉴权探活行为。
	HTTP client.Options
}

// ConfirmResult 是 Confirm 的业务响应：释放失败不是 RPC 错误，而是携带
// 本地副本完整性提示的 success=false，会话保留供下一次重试。
type ConfirmResult struct {
	// Message 是返回给 Python 客户端的提示文案。
	Message string
	// ReleaseFailed 表示 stop_experiment 未成功，会话仍保留在注册表中。
	ReleaseFailed bool
}

// Manager 维护 session_id 到会话的映射注册表，并统一编排 Run 级的生命周期。
//
// 注册表锁只用于定位、插入与移除会话；排空等待、stop_experiment 等可能阻塞的动作
// 一律在锁外执行。
type Manager struct {
	mu       sync.RWMutex
	sessions map[string]*Session
	released releasedCache
	factory  *client.Factory
	// healthReg 按站点共享探针，跨会话复用同一站点的健康状态。
	healthReg *health.Registry
}

// New 创建并初始化会话管理器实例。
func New(opts Options) *Manager {
	return &Manager{
		sessions:  make(map[string]*Session),
		released:  newReleasedCache(),
		factory:   client.NewFactory(opts.HTTP),
		healthReg: health.NewRegistry(nil),
	}
}

// install 生成 session_id 并把已装配完成的会话插入注册表。
// 插入前再次校验容量，避免装配期间其它会话占满槽位。
func (m *Manager) install(s *Session) (string, error) {
	id, err := uuid.NewV7()
	if err != nil {
		return "", fmt.Errorf("generate session id: %w", err)
	}
	sessionID := id.String()

	m.mu.Lock()
	defer m.mu.Unlock()
	if len(m.sessions) >= maxSessions {
		return "", fmt.Errorf("%w: limit %d", ErrTooManySessions, maxSessions)
	}
	m.sessions[sessionID] = s
	return sessionID, nil
}

// hasCapacity 在装配前做一次快速容量预检，避免为注定失败的会话申请外部资源。
func (m *Manager) hasCapacity() bool {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return len(m.sessions) < maxSessions
}

// DeliverRunFinish 执行指定会话的收尾序列（RUNNING → DRAINING）。
//
// 首次投递完成收尾；重复投递会从上次失败点续跑（FinishRecord 落盘、崩溃
// error_log 入队与 Intake.Close 都允许重试），因此失败后调用方可以安全重试。
// 收尾失败时返回原始失败信息，会话保留供 Confirm 续跑。
// 注册表锁只用于定位会话，刷盘、停 watcher 与关闭通道都在锁外执行。
func (m *Manager) DeliverRunFinish(ctx context.Context, sessionID string, rec *runv1.FinishRecord) error {
	if sessionID == "" {
		return ErrEmptySessionID
	}
	if rec == nil {
		return ErrInvalidFinishRecord
	}
	s, err := m.session(sessionID)
	if err != nil {
		return err
	}
	return s.finish(ctx, rec)
}

// State 获取指定会话的当前生命周期状态快照。
func (m *Manager) State(sessionID string) (RunState, error) {
	s, err := m.session(sessionID)
	if err != nil {
		return 0, err
	}
	return s.stateSnapshot(), nil
}

// Stats 返回指定会话的上传进度快照；计数由 tracker 维护，不持注册表锁。
func (m *Manager) Stats(sessionID string) (*operationv1.OperationStats, error) {
	s, err := m.session(sessionID)
	if err != nil {
		return nil, err
	}
	return s.snapshot(), nil
}

// confirmOKMessage 是成功释放时返回给 Python 客户端的固定文案。
const confirmOKMessage = "OK"

// Confirm 等待会话排空后释放其资源。
//
// 注册表只在定位会话和最终移除时持锁；排空等待与释放钩子都在锁外执行，
// 避免网络请求阻塞其它会话。并发 Confirm 由 per-session 的 confirmMu 串行化，
// 只有第一个调用会执行释放钩子，其余调用命中 released 缓存直接成功。
//
// 会话处于 FAILED 时先用新的 ctx 续跑被中断的收尾序列（FinishRecord 落盘、
// 崩溃 error_log 入队与 Intake.Close）；续跑仍失败就返回错误并保留会话，
// 不能在没有终态记录的情况下上报 stop_experiment。
//
// stop_experiment 失败不返回 RPC 错误，而是返回 ReleaseFailed=true 与按本地
// 副本完整性生成的提示文案，会话保留在注册表中供下一次 Confirm 重试。
func (m *Manager) Confirm(ctx context.Context, sessionID string) (ConfirmResult, error) {
	if sessionID == "" {
		return ConfirmResult{}, ErrEmptySessionID
	}
	s, released, err := m.lookup(sessionID)
	if err != nil {
		return ConfirmResult{}, err
	}
	if released {
		return ConfirmResult{Message: confirmOKMessage}, nil
	}
	if err = s.waitDrained(ctx); err != nil {
		return ConfirmResult{}, err
	}

	s.confirmMu.Lock()
	defer s.confirmMu.Unlock()

	// 等排空期间其它 Confirm 可能已经完成释放。
	s, released, err = m.lookup(sessionID)
	if err != nil {
		return ConfirmResult{}, err
	}
	if released {
		return ConfirmResult{Message: confirmOKMessage}, nil
	}
	switch s.stateSnapshot() {
	case RunStateFailed:
		// 收尾序列被中断：先续跑 FinishRecord 落盘、error_log 入队与 Intake.Close，
		// 全部成功后才允许上报终态并释放会话。
		if err = s.resumeFinish(ctx); err != nil {
			return ConfirmResult{}, err
		}
	case RunStateDrained:
		// Shutdown 只做上传排空、不写 FinishRecord，这类会话不能按成功路径释放。
		if !s.finishCommittedNow() {
			return ConfirmResult{}, fmt.Errorf("%w: run drained without a finish record", ErrRunFailed)
		}
	default:
		return ConfirmResult{}, ErrRunNotFinished
	}

	// 防御性检查：预算暂停会把记录留在本地而不进入上传管线。回填闭环（T4）
	// 完成前，任何暂停逗留都应阻止上报成功终态。
	if s.intake != nil {
		if paused, seq := s.intake.PauseState(); paused {
			return ConfirmResult{}, fmt.Errorf(
				"%w: uploads paused at store seq %d; records were not delivered", ErrRunFailed, seq)
		}
	}

	// 与 Python 对齐：排空成功后先停心跳，再上报实验终态；失败时心跳保持停止，
	// 会话保留供重试，避免重试期间继续发送心跳。
	if s.heartbeat != nil {
		s.heartbeat.stop(ctx)
		// 心跳等待可能耗尽 RPC 预算，ctx 已失效时不再发 stop_experiment。
		if err = ctx.Err(); err != nil {
			return ConfirmResult{}, err
		}
	}
	if s.releaseFn != nil {
		if err = s.releaseFn(ctx); err != nil {
			console.Warningf("confirm: stop experiment %q: %v", sessionID, err)
			return ConfirmResult{Message: s.confirmMessage(), ReleaseFailed: true}, nil
		}
	}
	if s.uploader != nil {
		_ = s.uploader.Close()
	}
	if s.cancel != nil {
		s.cancel()
	}
	s.closeHealth()

	m.mu.Lock()
	delete(m.sessions, sessionID)
	m.released.add(sessionID)
	m.mu.Unlock()
	return ConfirmResult{Message: confirmOKMessage}, nil
}

// Shutdown 有界回收全部活跃会话，供进程退出钩子调用。
//
// 每个会话并发执行「DRAINING → 停 watcher/heartbeat → 关本地 store →
// Intake.Close → RequestFinish → 等待排空」，不调用 stop_experiment；
// 总等待受 ctx 限制，超时后强制关闭所有 uploader 并告警。
func (m *Manager) Shutdown(ctx context.Context) error {
	m.mu.RLock()
	sessions := make([]*Session, 0, len(m.sessions))
	for _, s := range m.sessions {
		sessions = append(sessions, s)
	}
	m.mu.RUnlock()
	if len(sessions) == 0 {
		return nil
	}

	var wg sync.WaitGroup
	for _, s := range sessions {
		wg.Add(1)
		go func(s *Session) {
			defer wg.Done()
			s.shutdown(ctx)
		}(s)
	}
	done := make(chan struct{})
	go func() {
		wg.Wait()
		close(done)
	}()

	select {
	case <-done:
	case <-ctx.Done():
		console.Warningf("shutdown: %d active session(s) still draining, forcing uploader close", len(sessions))
	}

	// 无论是否在预算内排空，都释放 uploader 的 worker 协程；已排空时 Close 为 no-op。
	for _, s := range sessions {
		if s.uploader != nil {
			if err := s.uploader.Close(); err != nil {
				console.Warningf("shutdown: close uploader: %v", err)
			}
		}
	}

	select {
	case <-done:
		return nil
	case <-ctx.Done():
		return ctx.Err()
	}
}

// Close 释放 Manager 持有的进程级资源（HTTP 连接池），可重复调用。
func (m *Manager) Close() error {
	return m.factory.Close()
}

// lookup 在同一把读锁内区分「已释放」与「已注册」，避免两次读取之间被其它
// Confirm 移除而产生误报的 ErrSessionNotFound。
func (m *Manager) lookup(sessionID string) (*Session, bool, error) {
	m.mu.RLock()
	defer m.mu.RUnlock()
	if m.released.has(sessionID) {
		return nil, true, nil
	}
	s, ok := m.sessions[sessionID]
	if !ok {
		return nil, false, ErrSessionNotFound
	}
	return s, false, nil
}

func (m *Manager) session(sessionID string) (*Session, error) {
	if sessionID == "" {
		return nil, ErrEmptySessionID
	}
	m.mu.RLock()
	defer m.mu.RUnlock()
	s, ok := m.sessions[sessionID]
	if !ok {
		return nil, ErrSessionNotFound
	}
	return s, nil
}

func (m *Manager) isReleased(sessionID string) bool {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.released.has(sessionID)
}
