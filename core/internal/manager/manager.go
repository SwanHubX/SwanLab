package manager

import (
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
	// ErrRunNotFinished 表示会话尚未标记结束，不能提前确认释放。
	ErrRunNotFinished = errors.New("run finish has not been delivered for this handle")
)

// Manager 维护 session_id 到会话的映射注册表，并统一编排 Run 级的生命周期。
//
// 实现方式：
// 采用读写锁保护 sessions 映射表；所有针对会话的操作与释放均保证并发安全与原子性，
// 杜绝误释放未完成会话或产生竞态问题。
type Manager struct {
	mu       sync.RWMutex
	sessions map[string]*Session
}

// New 创建并初始化会话管理器实例。
func New() *Manager {
	return &Manager{sessions: make(map[string]*Session)}
}

// Start 创建新的 Run 会话，并返回 UUIDv7 形式的会话标识（session_id）。
// 每次调用都创建独立会话，不按实验身份去重：同一个 run_id 可能由多个进程上报
// （多 rank、resume 复用实验），它们的缓冲与上传状态必须互不影响。
func (m *Manager) Start() (string, error) {
	id, err := uuid.NewV7()
	if err != nil {
		return "", fmt.Errorf("generate session id: %w", err)
	}
	sessionID := id.String()

	m.mu.Lock()
	defer m.mu.Unlock()
	m.sessions[sessionID] = newSession()
	return sessionID, nil
}

// Finish 幂等标记指定会话为结束状态。
// 获取注册表读锁定位会话并在会话锁内更新状态，保证状态变更期间会话不被并发 Confirm 移除。
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
	s.finish()
	return nil
}

// State 获取指定会话的当前生命周期状态快照。
func (m *Manager) State(sessionID string) (RunState, error) {
	if sessionID == "" {
		return 0, ErrEmptySessionID
	}
	m.mu.RLock()
	defer m.mu.RUnlock()
	s, ok := m.sessions[sessionID]
	if !ok {
		return 0, ErrSessionNotFound
	}
	return s.stateSnapshot(), nil
}

// Confirm 原子校验会话是否已完成并释放会话资源。
// 若会话尚未完成，则拒绝释放并保留会话，供调用方后续重试。
func (m *Manager) Confirm(sessionID string) error {
	if sessionID == "" {
		return ErrEmptySessionID
	}
	m.mu.Lock()
	defer m.mu.Unlock()
	s, ok := m.sessions[sessionID]
	if !ok {
		return ErrSessionNotFound
	}
	if !s.finished() {
		return ErrRunNotFinished
	}
	delete(m.sessions, sessionID)
	return nil
}
