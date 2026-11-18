package syncer

import (
	"context"
	"errors"
	"fmt"
	"sync"

	"github.com/google/uuid"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

var (
	// ErrEmptySessionID 表示未提供有效的 sync 会话标识。
	ErrEmptySessionID = errors.New("sync session_id must not be empty")
	// ErrSessionNotFound 表示请求的 sync 会话不存在或已释放。
	ErrSessionNotFound = errors.New("unknown or released sync session")
	// ErrInvalidState 表示当前会话状态不允许执行请求的阶段操作。
	ErrInvalidState = errors.New("invalid sync session state")
	// ErrTooManySessions 表示活跃 sync 会话数量已达上限。
	ErrTooManySessions = errors.New("sync session registry is full")
	// ErrInvalidSettings 表示 sync 启动配置不合法。
	ErrInvalidSettings = errors.New("invalid sync settings")
	// ErrInvalidStartRecord 表示本地 run 文件缺失合法的 StartRecord。
	ErrInvalidStartRecord = errors.New("invalid sync start record")
)

const maxSyncSessions = 128

// Options 是 syncer.Manager 的构造参数。
type Options struct {
	HTTP client.Options
}

// StartOptions 描述 DeliverSyncStart 的输入参数。
type StartOptions struct {
	Settings  config.Settings
	Workspace string
	Project   string
	RunID     string
	Range     *corev1.SyncRange
}

// FlushResult 描述 DeliverSyncFlush 成功后返回的实验路径。
type FlushResult struct {
	Path string
}

// ConfirmResult 描述 ConfirmSyncFinish 的业务响应。
type ConfirmResult struct {
	Success bool
	Message string
}

// Manager 维护活跃的 sync 会话注册表并编排离线文件同步生命周期。
type Manager struct {
	mu       sync.RWMutex
	sessions map[string]*syncSession
	factory  *client.Factory
}

// NewManager 创建并初始化 sync 会话管理器。
func NewManager(opts Options) *Manager {
	return &Manager{
		sessions: make(map[string]*syncSession),
		factory:  client.NewFactory(opts.HTTP),
	}
}

// DeliverSyncStart 校验本地配置与范围参数、扫描 run-*.swanlab 首条 StartRecord 并创建 sync 会话。
// 本阶段不访问任何云端接口。
func (m *Manager) DeliverSyncStart(_ context.Context, opts StartOptions) (string, error) {
	if !m.hasCapacity() {
		return "", fmt.Errorf("%w: limit %d", ErrTooManySessions, maxSyncSessions)
	}
	filter, err := FromProto(opts.Range)
	if err != nil {
		return "", err
	}
	s, err := newSyncSession(opts, filter)
	if err != nil {
		return "", err
	}
	id, err := m.install(s)
	if err != nil {
		return "", err
	}
	return id, nil
}

// DeliverSyncFlush 完成云端实验创建/恢复与摘要水位初始化，并启动本地文件读取与上传管线。
func (m *Manager) DeliverSyncFlush(ctx context.Context, sessionID string) (FlushResult, error) {
	s, err := m.session(sessionID)
	if err != nil {
		return FlushResult{}, err
	}
	res, err := s.flush(ctx, m.factory)
	if err != nil {
		return FlushResult{}, err
	}
	return res, nil
}

// Stats 返回指定 sync 会话的上传进度快照。
func (m *Manager) Stats(sessionID string) (*operationv1.OperationStats, error) {
	s, err := m.session(sessionID)
	if err != nil {
		return nil, err
	}
	return s.stats(), nil
}

// ConfirmSyncFinish 等待本地扫描与上传全部排空，并按本地 FinishRecord 与新旧实验语义收尾。
func (m *Manager) ConfirmSyncFinish(ctx context.Context, sessionID string) (ConfirmResult, error) {
	s, err := m.session(sessionID)
	if err != nil {
		return ConfirmResult{}, err
	}
	res, err := s.confirm(ctx)
	if err != nil {
		return ConfirmResult{}, err
	}
	m.remove(sessionID)
	return res, nil
}

// Disconnect 在客户端连接断开后终止对应 sync 会话：取消扫描与上传但不上报云端终态。
func (m *Manager) Disconnect(ctx context.Context, sessionID string) error {
	if sessionID == "" {
		return ErrEmptySessionID
	}
	s := m.take(sessionID)
	if s == nil {
		return nil
	}
	return s.disconnect(ctx)
}

// Shutdown 在服务关闭时中止所有活跃 sync 会话（不向云端上报终态）。
func (m *Manager) Shutdown(ctx context.Context) error {
	m.mu.Lock()
	active := make([]*syncSession, 0, len(m.sessions))
	for id, s := range m.sessions {
		active = append(active, s)
		delete(m.sessions, id)
	}
	m.mu.Unlock()

	var errs []error
	for _, s := range active {
		if err := s.disconnect(ctx); err != nil {
			errs = append(errs, err)
		}
	}
	return errors.Join(errs...)
}

// Close 释放管理器持有的会话与传输资源。
func (m *Manager) Close() error {
	err := m.Shutdown(context.Background())
	if err != nil {
		console.Warningf("syncer shutdown: %v", err)
	}
	return err
}

func (m *Manager) hasCapacity() bool {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return len(m.sessions) < maxSyncSessions
}

func (m *Manager) install(s *syncSession) (string, error) {
	rawID, err := uuid.NewV7()
	if err != nil {
		return "", fmt.Errorf("generate sync session id: %w", err)
	}
	id := rawID.String()

	m.mu.Lock()
	defer m.mu.Unlock()
	if len(m.sessions) >= maxSyncSessions {
		return "", fmt.Errorf("%w: limit %d", ErrTooManySessions, maxSyncSessions)
	}
	m.sessions[id] = s
	return id, nil
}

func (m *Manager) session(id string) (*syncSession, error) {
	if id == "" {
		return nil, ErrEmptySessionID
	}
	m.mu.RLock()
	s, ok := m.sessions[id]
	m.mu.RUnlock()
	if !ok {
		return nil, fmt.Errorf("%w: %q", ErrSessionNotFound, id)
	}
	return s, nil
}

func (m *Manager) take(id string) *syncSession {
	m.mu.Lock()
	defer m.mu.Unlock()
	s := m.sessions[id]
	delete(m.sessions, id)
	return s
}

func (m *Manager) remove(id string) {
	m.mu.Lock()
	delete(m.sessions, id)
	m.mu.Unlock()
}
