package manager

import (
	"crypto/rand"
	"encoding/base64"
	"errors"
	"fmt"
	"sync"
)

var (
	// ErrEmptyHandle 表示未提供有效的 run_handle。
	ErrEmptyHandle = errors.New("run_handle must not be empty")
	// ErrRunNotFound 表示请求的 run_handle 不存在或已被确认释放。
	ErrRunNotFound = errors.New("unknown or released run_handle")
	// ErrRunNotFinished 表示会话尚未标记结束，不能提前确认释放。
	ErrRunNotFinished = errors.New("run finish has not been delivered for this handle")
)

// Manager 维护 run_handle 到会话的映射注册表，并统一编排 Run 级的生命周期。
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

// Start 创建新的 Run 会话，并返回随机生成的安全路由句柄（run_handle）。
func (m *Manager) Start() (string, error) {
	buf := make([]byte, 32)
	if _, err := rand.Read(buf); err != nil {
		return "", fmt.Errorf("generate run handle: %w", err)
	}
	handle := base64.RawURLEncoding.EncodeToString(buf)

	m.mu.Lock()
	defer m.mu.Unlock()
	if _, exists := m.sessions[handle]; exists {
		return "", errors.New("run handle collision")
	}
	m.sessions[handle] = newSession()
	return handle, nil
}

// Finish 幂等标记指定会话为结束状态。
// 获取注册表读锁定位会话并在会话锁内更新状态，保证状态变更期间会话不被并发 Confirm 移除。
func (m *Manager) Finish(handle string) error {
	if handle == "" {
		return ErrEmptyHandle
	}
	m.mu.RLock()
	defer m.mu.RUnlock()
	s, ok := m.sessions[handle]
	if !ok {
		return ErrRunNotFound
	}
	s.finish()
	return nil
}

// State 获取指定会话的当前生命周期状态快照。
func (m *Manager) State(handle string) (RunState, error) {
	if handle == "" {
		return 0, ErrEmptyHandle
	}
	m.mu.RLock()
	defer m.mu.RUnlock()
	s, ok := m.sessions[handle]
	if !ok {
		return 0, ErrRunNotFound
	}
	return s.stateSnapshot(), nil
}

// Confirm 原子校验会话是否已完成并释放会话资源。
// 若会话尚未完成，则拒绝释放并保留会话，供调用方后续重试。
func (m *Manager) Confirm(handle string) error {
	if handle == "" {
		return ErrEmptyHandle
	}
	m.mu.Lock()
	defer m.mu.Unlock()
	s, ok := m.sessions[handle]
	if !ok {
		return ErrRunNotFound
	}
	if !s.finished() {
		return ErrRunNotFinished
	}
	delete(m.sessions, handle)
	return nil
}
