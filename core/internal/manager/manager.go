package manager

import (
	"context"
	"errors"
	"fmt"
	"sync"

	"github.com/google/uuid"
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
	// ErrTooManySessions 表示会话数已达上限，需要先释放已结束的会话。
	ErrTooManySessions = errors.New("session registry is full")
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

// Manager 维护 session_id 到会话的映射注册表，并统一编排 Run 级的生命周期。
//
// 注册表锁只用于定位、插入与移除会话；排空等待、stop_experiment 等可能阻塞的动作
// 一律在锁外执行。
type Manager struct {
	mu       sync.RWMutex
	sessions map[string]*Session
	released releasedCache
}

// New 创建并初始化会话管理器实例。
func New() *Manager {
	return &Manager{
		sessions: make(map[string]*Session),
		released: newReleasedCache(),
	}
}

// Start 创建新的 Run 会话，并返回 UUIDv7 形式的会话标识（session_id）。
// 每次调用都创建独立会话，不按实验身份去重：同一个 run_id 可能由多个进程上报
// （多 rank、resume 复用实验），它们的缓冲与上传状态必须互不影响。
func (m *Manager) Start() (string, error) {
	return m.start(sessionOptions{})
}

// start 是 Start 的内部实现，允许测试注入排空与释放钩子。
func (m *Manager) start(opts sessionOptions) (string, error) {
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
	m.sessions[sessionID] = newSession(opts)
	return sessionID, nil
}

// Finish 幂等标记指定会话进入排空阶段（RUNNING → DRAINING）。
// 已排空或已在排空时返回成功；收尾失败时返回原始失败信息。
// 持注册表读锁执行切换，保证不会对已被并发 Confirm 释放的会话返回成功。
func (m *Manager) Finish(sessionID string) error {
	if sessionID == "" {
		return ErrEmptySessionID
	}
	m.mu.RLock()
	defer m.mu.RUnlock()
	s, ok := m.sessions[sessionID]
	if !ok {
		return ErrSessionNotFound
	}
	return s.beginDrain()
}

// State 获取指定会话的当前生命周期状态快照。
func (m *Manager) State(sessionID string) (RunState, error) {
	s, err := m.session(sessionID)
	if err != nil {
		return 0, err
	}
	return s.stateSnapshot(), nil
}

// Confirm 等待会话排空后释放其资源。
//
// 注册表只在定位会话和最终移除时持锁；排空等待与释放钩子都在锁外执行，
// 避免网络请求阻塞其它会话。并发 Confirm 由 per-session 的 confirmMu 串行化，
// 只有第一个调用会执行释放钩子，其余调用命中 released 缓存直接成功。
func (m *Manager) Confirm(ctx context.Context, sessionID string) error {
	if sessionID == "" {
		return ErrEmptySessionID
	}
	s, released, err := m.lookup(sessionID)
	if err != nil {
		return err
	}
	if released {
		return nil
	}
	if err := s.waitDrained(ctx); err != nil {
		return err
	}

	s.confirmMu.Lock()
	defer s.confirmMu.Unlock()

	// 等排空期间其它 Confirm 可能已经完成释放。
	s, released, err = m.lookup(sessionID)
	if err != nil {
		return err
	}
	if released {
		return nil
	}
	switch s.stateSnapshot() {
	case RunStateFailed:
		return s.failure()
	case RunStateDrained:
	default:
		return ErrRunNotFinished
	}

	if s.releaseFn != nil {
		if err := s.releaseFn(ctx); err != nil {
			return err
		}
	}

	m.mu.Lock()
	delete(m.sessions, sessionID)
	m.released.add(sessionID)
	m.mu.Unlock()
	return nil
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
