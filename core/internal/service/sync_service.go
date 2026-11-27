package service

import (
	"context"
	"errors"
	"sync"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"
	"google.golang.org/protobuf/types/known/emptypb"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/server"
	"github.com/swanhubx/swanlab/core/internal/syncer"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

const syncOKMessage = "success"

// SyncManager 是 SyncService 依赖的 sync 会话编排接口，由 syncer.Manager 实现；
// gRPC 单测可注入 Fake，避免每个 RPC 用例都依赖磁盘扫描与云端交互。
type SyncManager interface {
	DeliverSyncStart(context.Context, syncer.StartOptions) (string, error)
	DeliverSyncFlush(context.Context, string) (syncer.FlushResult, error)
	Stats(string) (*operationv1.OperationStats, error)
	ConfirmSyncFinish(context.Context, string) (syncer.ConfirmResult, error)
	Disconnect(context.Context, string) error
}

var _ SyncManager = (*syncer.Manager)(nil)

// SyncService 实现 CoreSyncServiceServer 契约，按 gRPC 连接绑定并驱动本地 Run 离线同步会话。
type SyncService struct {
	corev1.UnimplementedCoreSyncServiceServer

	lifecycle *server.Lifecycle
	syncs     SyncManager

	connMu      sync.Mutex
	pendingConn map[uint64]struct{}
	byConn      map[uint64]string
	cleanups    *disconnectPool
}

// NewSyncService 创建 CoreSyncService gRPC 接入实例，关联服务生命周期与 syncer 管理器。
func NewSyncService(controller *server.Controller, syncManager SyncManager) *SyncService {
	var lifecycle *server.Lifecycle
	if controller != nil {
		lifecycle = controller.Lifecycle()
	}
	return &SyncService{
		lifecycle:   lifecycle,
		syncs:       syncManager,
		pendingConn: make(map[uint64]struct{}),
		byConn:      make(map[uint64]string),
		cleanups:    newDisconnectPool(defaultDisconnectGrace),
	}
}

// Register 将 CoreSyncService 注册到指定的 gRPC Server。
func (s *SyncService) Register(g *grpc.Server) {
	corev1.RegisterCoreSyncServiceServer(g, s)
}

func (s *SyncService) requireActive() error {
	if s.lifecycle == nil {
		return nil
	}
	switch s.lifecycle.Get() {
	case server.StateStopping, server.StateClosed:
		return status.Error(codes.FailedPrecondition, "service is stopping or closed")
	default:
		return nil
	}
}

// DeliverSyncStart 校验启动配置与本地 run-*.swanlab 文件头，为当前连接创建并绑定同步会话。
func (s *SyncService) DeliverSyncStart(
	ctx context.Context,
	req *corev1.DeliverSyncStartRequest,
) (*corev1.DeliverSyncStartResponse, error) {
	if err := s.requireActive(); err != nil {
		return nil, err
	}
	if req == nil || req.GetCoreSettings() == nil {
		return nil, status.Error(codes.InvalidArgument, "core_settings must not be nil")
	}
	settings, err := protoToSettings(req.GetCoreSettings())
	if err != nil {
		return nil, status.Error(codes.InvalidArgument, err.Error())
	}

	connID := syncConnID(ctx)
	if resErr := s.reserveConn(ctx, connID); resErr != nil {
		return nil, resErr
	}

	sid, err := s.syncs.DeliverSyncStart(ctx, syncer.StartOptions{
		Settings:  settings,
		Workspace: req.GetWorkspace(),
		Project:   req.GetProject(),
		RunID:     req.GetId(),
		Range:     req.GetRange(),
	})
	if err != nil {
		s.releasePending(connID)
		if isSyncStartRejection(err) {
			return &corev1.DeliverSyncStartResponse{Success: false, Message: err.Error()}, nil
		}
		return nil, syncStatus(err)
	}

	if !s.commitConnSession(ctx, connID, sid) {
		s.cleanupDisconnected(sid)
		return nil, status.Error(codes.Canceled, "client connection closed during DeliverSyncStart")
	}
	return &corev1.DeliverSyncStartResponse{Success: true, Message: syncOKMessage}, nil
}

// DeliverSyncFlush 初始化云端实验与摘要水位，并启动当前连接对应 sync 会话的后台扫描与上传。
func (s *SyncService) DeliverSyncFlush(
	ctx context.Context,
	_ *emptypb.Empty,
) (*corev1.DeliverSyncFlushResponse, error) {
	if err := s.requireActive(); err != nil {
		return nil, err
	}
	sid, ok := s.boundSession(syncConnID(ctx))
	if !ok {
		return nil, status.Error(codes.FailedPrecondition, "sync session not started; call DeliverSyncStart first")
	}

	res, err := s.syncs.DeliverSyncFlush(ctx, sid)
	if err != nil {
		if errors.Is(err, syncer.ErrSessionNotFound) || errors.Is(err, syncer.ErrInvalidState) {
			return nil, status.Error(codes.FailedPrecondition, err.Error())
		}
		if ctxErr := syncContextStatus(ctx, err); ctxErr != nil {
			return nil, ctxErr
		}
		return &corev1.DeliverSyncFlushResponse{Success: false, Message: err.Error()}, nil
	}
	return &corev1.DeliverSyncFlushResponse{
		Success: true,
		Message: syncOKMessage,
		Path:    res.Path,
	}, nil
}

// GetOperationStats 返回当前连接绑定 sync 会话的实时上传与运行状态快照。
func (s *SyncService) GetOperationStats(
	ctx context.Context,
	_ *emptypb.Empty,
) (*corev1.GetOperationStatsResponse, error) {
	if err := s.requireActive(); err != nil {
		return nil, err
	}
	sid, ok := s.boundSession(syncConnID(ctx))
	if !ok {
		return &corev1.GetOperationStatsResponse{
			Success: false,
			Message: "sync session not started",
		}, nil
	}
	stats, err := s.syncs.Stats(sid)
	if err != nil {
		if errors.Is(err, syncer.ErrSessionNotFound) {
			return &corev1.GetOperationStatsResponse{Success: false, Message: err.Error()}, nil
		}
		return nil, syncStatus(err)
	}
	return &corev1.GetOperationStatsResponse{
		Success: true,
		Stats:   stats,
	}, nil
}

// ConfirmSyncFinish 等待扫描与上传排空完成并按运行终态收尾，然后释放当前连接的 sync 会话绑定。
func (s *SyncService) ConfirmSyncFinish(
	ctx context.Context,
	_ *emptypb.Empty,
) (*corev1.ConfirmSyncFinishResponse, error) {
	if err := s.requireActive(); err != nil {
		return nil, err
	}
	connID := syncConnID(ctx)
	sid, ok := s.boundSession(connID)
	if !ok {
		return nil, status.Error(codes.FailedPrecondition, "sync session not started; call DeliverSyncStart first")
	}

	res, err := s.syncs.ConfirmSyncFinish(ctx, sid)
	if err != nil {
		if errors.Is(err, syncer.ErrInvalidState) {
			return nil, status.Error(codes.FailedPrecondition, err.Error())
		}
		if ctxErr := syncContextStatus(ctx, err); ctxErr != nil {
			return nil, ctxErr
		}
		return nil, syncStatus(err)
	}
	s.unbindConnSession(connID, sid)
	msg := res.Message
	if msg == "" {
		msg = syncOKMessage
	}
	return &corev1.ConfirmSyncFinishResponse{
		Success: res.Success,
		Message: msg,
	}, nil
}

// HandleConnDisconnect 在传输连接断开后异步取消并释放该连接绑定但尚未完成确认的 sync 会话，
// 绝不在未完成 ConfirmSyncFinish 时向云端误报 StopExperiment。
func (s *SyncService) HandleConnDisconnect(connID uint64) {
	if connID == 0 {
		return
	}
	s.connMu.Lock()
	delete(s.pendingConn, connID)
	sid, ok := s.byConn[connID]
	if ok {
		delete(s.byConn, connID)
	}
	s.connMu.Unlock()
	if !ok {
		return
	}
	s.cleanupDisconnected(sid)
}

// WaitDisconnectCleanups 封闭新的断连清理派发并等待在途清理协程结束；超时后取消清理上下文。
func (s *SyncService) WaitDisconnectCleanups(ctx context.Context) error {
	return s.cleanups.wait(ctx)
}

func syncConnID(ctx context.Context) uint64 {
	connID, ok := server.ConnIDFromContext(ctx)
	if !ok {
		return 0
	}
	return connID
}

func (s *SyncService) reserveConn(ctx context.Context, connID uint64) error {
	s.connMu.Lock()
	defer s.connMu.Unlock()
	if s.cleanups.isClosed() || server.IsConnClosed(ctx) {
		return status.Error(codes.Canceled, "client connection closed")
	}
	if _, exists := s.byConn[connID]; exists {
		return status.Error(codes.FailedPrecondition, "sync session already active on this connection")
	}
	if _, pending := s.pendingConn[connID]; pending {
		return status.Error(codes.FailedPrecondition, "sync session start already in progress on this connection")
	}
	s.pendingConn[connID] = struct{}{}
	return nil
}

func (s *SyncService) releasePending(connID uint64) {
	s.connMu.Lock()
	delete(s.pendingConn, connID)
	s.connMu.Unlock()
}

func (s *SyncService) commitConnSession(ctx context.Context, connID uint64, sessionID string) bool {
	if sessionID == "" {
		return true
	}
	s.connMu.Lock()
	defer s.connMu.Unlock()
	_, stillPending := s.pendingConn[connID]
	delete(s.pendingConn, connID)
	if !stillPending || s.cleanups.isClosed() || server.IsConnClosed(ctx) {
		return false
	}
	s.byConn[connID] = sessionID
	return true
}

func (s *SyncService) boundSession(connID uint64) (string, bool) {
	s.connMu.Lock()
	defer s.connMu.Unlock()
	sid, ok := s.byConn[connID]
	return sid, ok && sid != ""
}

func (s *SyncService) unbindConnSession(connID uint64, sessionID string) {
	s.connMu.Lock()
	defer s.connMu.Unlock()
	if cur, ok := s.byConn[connID]; ok && cur == sessionID {
		delete(s.byConn, connID)
	}
}

func (s *SyncService) cleanupDisconnected(sessionID string) {
	if sessionID == "" {
		return
	}
	s.cleanups.goRun(func(ctx context.Context) {
		if err := s.syncs.Disconnect(ctx, sessionID); err != nil {
			console.Warningf("sync disconnect cleanup for session %q failed: %v", sessionID, err)
		}
	})
}

func isSyncStartRejection(err error) bool {
	return errors.Is(err, syncer.ErrInvalidSettings) ||
		errors.Is(err, syncer.ErrInvalidStartRecord) ||
		errors.Is(err, syncer.ErrInvalidSyncRange)
}

func syncContextStatus(ctx context.Context, err error) error {
	if errors.Is(err, context.DeadlineExceeded) || errors.Is(ctx.Err(), context.DeadlineExceeded) {
		return status.Error(codes.DeadlineExceeded, err.Error())
	}
	if errors.Is(err, context.Canceled) || errors.Is(ctx.Err(), context.Canceled) {
		return status.Error(codes.Canceled, err.Error())
	}
	return nil
}

func syncStatus(err error) error {
	switch {
	case errors.Is(err, syncer.ErrEmptySessionID):
		return status.Error(codes.InvalidArgument, err.Error())
	case errors.Is(err, syncer.ErrSessionNotFound):
		return status.Error(codes.NotFound, err.Error())
	case errors.Is(err, syncer.ErrInvalidState):
		return status.Error(codes.FailedPrecondition, err.Error())
	case errors.Is(err, syncer.ErrInvalidSettings),
		errors.Is(err, syncer.ErrInvalidStartRecord),
		errors.Is(err, syncer.ErrInvalidSyncRange):
		return status.Error(codes.InvalidArgument, err.Error())
	case errors.Is(err, syncer.ErrTooManySessions):
		return status.Error(codes.ResourceExhausted, err.Error())
	case errors.Is(err, context.DeadlineExceeded):
		return status.Error(codes.DeadlineExceeded, err.Error())
	case errors.Is(err, context.Canceled):
		return status.Error(codes.Canceled, err.Error())
	default:
		return status.Error(codes.Internal, err.Error())
	}
}
