// Package service 提供 swanlab-core 的 gRPC 协议接入层实现。
//
// 职责：
// - 实现 gRPC 契约（CoreServiceServer），校验并解析 Protobuf 请求；
// - 拦截未就绪服务的业务调用，前置统一做状态检查；
// - 将内部领域错误（manager）映射为标准的 gRPC 状态码（InvalidArgument、NotFound 等）；
// - 将服务退出指令委派给 internal/server 控制器执行。
package service

import (
	"context"
	"errors"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"
	"google.golang.org/protobuf/types/known/emptypb"

	"github.com/swanhubx/swanlab/core/internal/manager"
	"github.com/swanhubx/swanlab/core/internal/server"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

// CoreService 实现生成的 CoreServiceServer 契约，作为 gRPC 请求的接入网关。
type CoreService struct {
	corev1.UnimplementedCoreServiceServer

	controller *server.Controller
	lifecycle  *server.Lifecycle
	manager    *manager.Manager
}

// NewCoreService 创建 gRPC 接入服务实例，关联服务生命周期控制器与 Run 会话管理器。
func NewCoreService(controller *server.Controller, runManager *manager.Manager) *CoreService {
	return &CoreService{
		controller: controller,
		lifecycle:  controller.Lifecycle(),
		manager:    runManager,
	}
}

// Register 将当前服务注册到指定的 gRPC Server。
func (s *CoreService) Register(g *grpc.Server) {
	corev1.RegisterCoreServiceServer(g, s)
}

// requireReady 检查服务是否处于就绪状态，未就绪时返回 FailedPrecondition 错误。
func (s *CoreService) requireReady() error {
	if s.lifecycle.Ready() {
		return nil
	}
	return status.Error(codes.FailedPrecondition, "service is not ready; call SpinupService first")
}

// rejectUpsert 统一拦截数据上传 RPC，在底层存储与传输通道接入前明确返回 Unimplemented。
func (s *CoreService) rejectUpsert() error {
	if err := s.requireReady(); err != nil {
		return err
	}
	return status.Error(codes.Unimplemented, "upsert data path is not implemented in this core build")
}

// SpinupService 执行服务启动握手，将服务状态置为就绪（Ready）。
func (s *CoreService) SpinupService(
	_ context.Context,
	_ *corev1.SpinupServiceRequest,
) (*corev1.SpinupServiceResponse, error) {
	if !s.lifecycle.Spinup() {
		return nil, status.Error(codes.FailedPrecondition, "service is stopping or closed")
	}
	return &corev1.SpinupServiceResponse{}, nil
}

// TeardownService 请求关闭整个 Core 服务，异步触发控制器的关闭流程。
func (s *CoreService) TeardownService(
	_ context.Context,
	_ *corev1.TeardownServiceRequest,
) (*corev1.TeardownServiceResponse, error) {
	s.controller.Shutdown("teardown")
	return &corev1.TeardownServiceResponse{}, nil
}

// DeliverRunStart 接收并校验启动记录，创建新的 Run 会话并返回对应的 run_handle。
func (s *CoreService) DeliverRunStart(
	_ context.Context,
	req *corev1.DeliverRunStartRequest,
) (*corev1.DeliverRunStartResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if req.GetStartRecord() == nil {
		return nil, status.Error(codes.InvalidArgument, "start_record must not be nil")
	}
	handle, err := s.manager.Start()
	if err != nil {
		return nil, status.Errorf(codes.Internal, "start run session: %v", err)
	}
	return &corev1.DeliverRunStartResponse{
		Success:       true,
		Run:           req.GetStartRecord(),
		NewExperiment: true,
		RunHandle:     handle,
	}, nil
}

func (s *CoreService) UpsertColumns(context.Context, *corev1.UpsertColumnsRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *CoreService) UpsertScalars(context.Context, *corev1.UpsertScalarsRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *CoreService) UpsertMedia(context.Context, *corev1.UpsertMediaRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *CoreService) UpsertLogs(context.Context, *corev1.UpsertLogsRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *CoreService) UpsertSaves(context.Context, *corev1.UpsertSavesRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

// DeliverRunFinish 标记指定 run_handle 的会话结束。
func (s *CoreService) DeliverRunFinish(
	_ context.Context,
	req *corev1.DeliverRunFinishRequest,
) (*corev1.DeliverRunFinishResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if err := s.manager.Finish(req.GetRunHandle()); err != nil {
		return nil, managerStatus(err)
	}
	return &corev1.DeliverRunFinishResponse{Success: true}, nil
}

// GetOperationStats 查询指定会话的当前运行与排空状态。
func (s *CoreService) GetOperationStats(
	_ context.Context,
	req *corev1.GetOperationStatsRequest,
) (*corev1.GetOperationStatsResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	runState, err := s.manager.State(req.GetRunHandle())
	if err != nil {
		return nil, managerStatus(err)
	}
	state := operationv1.CoreState_CORE_STATE_RUNNING
	if runState == manager.RunStateFinished {
		state = operationv1.CoreState_CORE_STATE_FINISHED
	}
	return &corev1.GetOperationStatsResponse{
		Success: true,
		Stats:   &operationv1.OperationStats{State: state},
	}, nil
}

// ConfirmRunFinish 确认会话已完成并释放相关资源。
func (s *CoreService) ConfirmRunFinish(
	_ context.Context,
	req *corev1.ConfirmRunFinishRequest,
) (*corev1.ConfirmRunFinishResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if err := s.manager.Confirm(req.GetRunHandle()); err != nil {
		return nil, managerStatus(err)
	}
	return &corev1.ConfirmRunFinishResponse{Success: true}, nil
}

// managerStatus 将 manager 层的领域错误转换为对应的 gRPC 状态码。
func managerStatus(err error) error {
	switch {
	case errors.Is(err, manager.ErrEmptyHandle):
		return status.Error(codes.InvalidArgument, err.Error())
	case errors.Is(err, manager.ErrRunNotFound):
		return status.Error(codes.NotFound, err.Error())
	case errors.Is(err, manager.ErrRunNotFinished):
		return status.Error(codes.FailedPrecondition, err.Error())
	default:
		return status.Error(codes.Internal, err.Error())
	}
}
