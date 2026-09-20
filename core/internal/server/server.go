// Package server 提供 swanlab-core 的 gRPC 服务端实现。
//
// 生命周期（P0 骨架）：
//   - auth interceptor 校验每个 RPC 的 auth token（port-file 回报）；
//   - SpinupService 是 owner-only、幂等的 READY 屏障；
//   - 五类 Upsert* 在 store/transport 实现前统一返回 UNIMPLEMENTED，不假成功；
//   - DeliverRunStart 创建 run_handle 会话，finish/stats/confirm 按 handle 路由；
//   - ConfirmRunFinish 只释放对应会话，绝不关闭 gRPC Server；
//   - TeardownService 校验 owner token，触发 controller 的统一关闭路径。
package server

import (
	"context"
	"crypto/hmac"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"
	"google.golang.org/protobuf/types/known/emptypb"

	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

// Service 实现 CoreService。
//
// 除 Teardown/Spinup 与 run 生命周期 RPC 外，五类业务 upsert 由嵌入的
// Unimplemented 接管并显式改写为 UNIMPLEMENTED（在 READY 门控之后）。
type Service struct {
	corev1.UnimplementedCoreServiceServer

	ownerToken string
	authToken  string
	controller *Controller
	lc         *Lifecycle
	registry   *sessionRegistry
}

// NewService 创建服务实例。ownerToken 为 spawn owner 通过私有文件传入的
// 服务级令牌，是唯一允许触发 Spinup/Teardown 的凭证；authToken 为
// port-file 回报给所有 client 的 RPC 鉴权令牌。
func NewService(ownerToken, authToken string, controller *Controller) *Service {
	return &Service{
		ownerToken: ownerToken,
		authToken:  authToken,
		controller: controller,
		lc:         controller.Lifecycle(),
		registry:   newSessionRegistry(),
	}
}

// Register 将服务注册到 gRPC Server。
func (s *Service) Register(g *grpc.Server) {
	corev1.RegisterCoreServiceServer(g, s)
}

// verifyOwner 常量时间校验 owner token；token 值不得出现在日志或错误消息中。
func (s *Service) verifyOwner(token string) error {
	if s.ownerToken == "" || !hmac.Equal([]byte(token), []byte(s.ownerToken)) {
		return status.Error(codes.PermissionDenied, "invalid owner token")
	}
	return nil
}

// requireReady 拒绝 READY 之前的 run 级 RPC。
func (s *Service) requireReady() error {
	if s.lc.Ready() {
		return nil
	}
	return status.Error(codes.FailedPrecondition, "service is not ready; call SpinupService first")
}

// rejectUpsert 是五类 upsert 的统一应答：READY 门控通过后必须返回
// UNIMPLEMENTED——store/transport 未实现前不得用 Empty 成功假接收。
func (s *Service) rejectUpsert() error {
	if err := s.requireReady(); err != nil {
		return err
	}
	return status.Error(codes.Unimplemented, "upsert data path is not implemented in this core build")
}

// SpinupService 完成 owner 校验并把服务置为 READY，幂等；
// 已进入 STOPPING/CLOSED 时拒绝。
func (s *Service) SpinupService(_ context.Context, req *corev1.SpinupServiceRequest) (*corev1.SpinupServiceResponse, error) {
	if err := s.verifyOwner(req.GetOwnerToken()); err != nil {
		return nil, err
	}
	if !s.lc.Spinup() {
		return nil, status.Error(codes.FailedPrecondition, "service is stopping or closed")
	}
	return &corev1.SpinupServiceResponse{}, nil
}

// TeardownService 关闭整个服务进程，仅接受正确的 owner token。
// 校验使用常量时间比较；token 值不得出现在日志或错误消息中。
func (s *Service) TeardownService(_ context.Context, req *corev1.TeardownServiceRequest) (*corev1.TeardownServiceResponse, error) {
	if err := s.verifyOwner(req.GetOwnerToken()); err != nil {
		return nil, err
	}
	// 异步触发统一关闭路径，保证本响应先于连接关闭送达调用方。
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
	// P0 骨架：仅回显最小有效结果，不落盘、不上传、不访问云端。
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

// DeliverRunFinish 只结束对应会话的 run，服务保持存活。
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

// ConfirmRunFinish 确认对应会话已排空并释放资源，绝不关闭 gRPC Server。
func (s *Service) ConfirmRunFinish(_ context.Context, req *corev1.ConfirmRunFinishRequest) (*corev1.ConfirmRunFinishResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if _, err := s.registry.lookup(req.GetRunHandle()); err != nil {
		return nil, err
	}
	s.registry.release(req.GetRunHandle())
	return &corev1.ConfirmRunFinishResponse{Success: true}, nil
}
