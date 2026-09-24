// Package server 提供 swanlab-core 的 gRPC 服务端实现。
//
// 生命周期（骨架）：
//   - 不使用应用层 token；信任边界是本地 transport 与文件系统权限
//     （UDS socket、owner-only runtime 目录、0600 port-file），能连接端点的
//     client 可调用全部 RPC，包括 Teardown；
//   - SpinupService 是幂等的 READY 屏障；
//   - 五类 Upsert* 在 store/transport 实现前返回 UNIMPLEMENTED；
//   - DeliverRunStart 创建 run_handle 会话，finish/stats/confirm 按 handle 路由；
//   - ConfirmRunFinish 在 finish 交付后原子释放对应会话，不关闭 gRPC Server；
//   - TeardownService 触发 controller 的关闭路径。
package server

import (
	"context"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"
	"google.golang.org/protobuf/types/known/emptypb"

	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

// Service 实现 CoreService。
//
// 除 Spinup/Teardown 与 run 生命周期 RPC 外，五类业务 upsert 由嵌入的
// Unimplemented 改写为 UNIMPLEMENTED（READY 门控之后）。
type Service struct {
	corev1.UnimplementedCoreServiceServer

	controller *Controller
	lc         *Lifecycle
	registry   *sessionRegistry
}

// NewService 创建服务实例。
func NewService(controller *Controller) *Service {
	return &Service{
		controller: controller,
		lc:         controller.Lifecycle(),
		registry:   newSessionRegistry(),
	}
}

// Register 将服务注册到 gRPC Server。
func (s *Service) Register(g *grpc.Server) {
	corev1.RegisterCoreServiceServer(g, s)
}

// requireReady 拒绝 READY 之前的 run 级 RPC。
func (s *Service) requireReady() error {
	if s.lc.Ready() {
		return nil
	}
	return status.Error(codes.FailedPrecondition, "service is not ready; call SpinupService first")
}

// rejectUpsert 是五类 upsert 的应答：READY 门控通过后返回
// UNIMPLEMENTED——store/transport 未实现前不得用 Empty 成功假接收。
func (s *Service) rejectUpsert() error {
	if err := s.requireReady(); err != nil {
		return err
	}
	return status.Error(codes.Unimplemented, "upsert data path is not implemented in this core build")
}

// SpinupService 把服务置为 READY，幂等；STOPPING/CLOSED 下拒绝。
func (s *Service) SpinupService(_ context.Context, _ *corev1.SpinupServiceRequest) (*corev1.SpinupServiceResponse, error) {
	if !s.lc.Spinup() {
		return nil, status.Error(codes.FailedPrecondition, "service is stopping or closed")
	}
	return &corev1.SpinupServiceResponse{}, nil
}

// TeardownService 关闭整个服务进程。无 owner 校验：能连接端点的 client
// 都可调用，SDK 约定由 owner 发起。异步触发关闭路径，响应先于连接关闭送达。
func (s *Service) TeardownService(_ context.Context, _ *corev1.TeardownServiceRequest) (*corev1.TeardownServiceResponse, error) {
	s.controller.Shutdown("teardown")
	return &corev1.TeardownServiceResponse{}, nil
}

// DeliverRunStart 创建 run 会话并返回非空 opaque run_handle。
func (s *Service) DeliverRunStart(_ context.Context, req *corev1.DeliverRunStartRequest) (*corev1.DeliverRunStartResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if req.GetStartRecord() == nil {
		return nil, status.Error(codes.InvalidArgument, "start_record must not be nil")
	}
	handle, _, err := s.registry.create()
	if err != nil {
		return nil, err
	}
	// 骨架：回显最小有效结果，不落盘、不上传、不访问云端。
	return &corev1.DeliverRunStartResponse{
		Success:       true,
		Run:           req.GetStartRecord(),
		NewExperiment: true,
		RunHandle:     handle,
	}, nil
}

func (s *Service) UpsertColumns(_ context.Context, _ *corev1.UpsertColumnsRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *Service) UpsertScalars(_ context.Context, _ *corev1.UpsertScalarsRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *Service) UpsertMedia(_ context.Context, _ *corev1.UpsertMediaRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *Service) UpsertLogs(_ context.Context, _ *corev1.UpsertLogsRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

func (s *Service) UpsertSaves(_ context.Context, _ *corev1.UpsertSavesRequest) (*emptypb.Empty, error) {
	return nil, s.rejectUpsert()
}

// DeliverRunFinish 结束对应会话的 run，服务保持存活。
func (s *Service) DeliverRunFinish(_ context.Context, req *corev1.DeliverRunFinishRequest) (*corev1.DeliverRunFinishResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	session, err := s.registry.lookup(req.GetRunHandle())
	if err != nil {
		return nil, err
	}
	session.finish()
	return &corev1.DeliverRunFinishResponse{Success: true}, nil
}

// GetOperationStats 返回对应会话的 CoreState 映射（run 级数据排空轴）。
func (s *Service) GetOperationStats(_ context.Context, req *corev1.GetOperationStatsRequest) (*corev1.GetOperationStatsResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	session, err := s.registry.lookup(req.GetRunHandle())
	if err != nil {
		return nil, err
	}
	return &corev1.GetOperationStatsResponse{
		Success: true,
		Stats:   &operationv1.OperationStats{State: session.coreState()},
	}, nil
}

// ConfirmRunFinish 确认对应会话已排空并原子释放资源，不关闭 gRPC Server。
// finish 未交付时返回 FAILED_PRECONDITION，会话保留。
func (s *Service) ConfirmRunFinish(_ context.Context, req *corev1.ConfirmRunFinishRequest) (*corev1.ConfirmRunFinishResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if err := s.registry.confirm(req.GetRunHandle()); err != nil {
		return nil, err
	}
	return &corev1.ConfirmRunFinishResponse{Success: true}, nil
}
