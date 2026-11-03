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
	"fmt"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"
	"google.golang.org/protobuf/types/known/emptypb"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager"
	"github.com/swanhubx/swanlab/core/internal/server"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	settingsv1 "github.com/swanhubx/swanlab/core/proto/swanlab/settings/core/v1"
)

// RunManager 是 service 依赖的会话编排接口，由 manager.Manager 实现；
// gRPC 单测可注入 Fake，避免每个 RPC 用例都依赖真实云端。
type RunManager interface {
	DeliverRunStart(context.Context, config.Settings, *runv1.StartRecord) (*manager.StartResult, error)
	Finish(string) error
	State(string) (manager.RunState, error)
	Confirm(context.Context, string) error
}

// CoreService 实现生成的 CoreServiceServer 契约，作为 gRPC 请求的接入网关。
type CoreService struct {
	corev1.UnimplementedCoreServiceServer

	controller *server.Controller
	lifecycle  *server.Lifecycle
	runs       RunManager
}

// NewCoreService 创建 gRPC 接入服务实例，关联服务生命周期控制器与 Run 会话管理器。
func NewCoreService(controller *server.Controller, runManager RunManager) *CoreService {
	return &CoreService{
		controller: controller,
		lifecycle:  controller.Lifecycle(),
		runs:       runManager,
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

// DeliverRunStart 接收并校验启动记录，创建新的 Run 会话并返回对应的 session_id。
// 每次调用创建一个独立会话：同一 run_id 可能来自多个进程（多 rank、resume 复用实验），
// 它们共享实验身份，但缓冲与上传状态必须各自独立。
// run_id 由 SDK 生成后透传，core 不补写。
func (s *CoreService) DeliverRunStart(
	ctx context.Context,
	req *corev1.DeliverRunStartRequest,
) (*corev1.DeliverRunStartResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	record := req.GetStartRecord()
	if record == nil {
		return nil, status.Error(codes.InvalidArgument, "start_record must not be nil")
	}
	settings, err := settingsFromProto(req.GetCoreSettings())
	if err != nil {
		return nil, status.Error(codes.InvalidArgument, err.Error())
	}
	result, err := s.runs.DeliverRunStart(ctx, settings, record)
	if err != nil {
		return nil, managerStatus(err)
	}
	if result.Path == "" {
		return nil, status.Error(codes.InvalidArgument, "project path is empty; online mode requires a path")
	}
	return &corev1.DeliverRunStartResponse{
		Success:          true,
		Run:              result.Run,
		Path:             result.Path,
		Name:             result.Name,
		GlobalStep:       result.GlobalStep,
		GlobalSystemStep: result.GlobalSystemStep,
		NewExperiment:    result.NewExperiment,
		SessionId:        result.SessionID,
	}, nil
}

// settingsFromProto 把 CoreSettings 一一映射为纯 Go 配置，并补齐 pydantic 未覆盖的下界。
func settingsFromProto(in *settingsv1.CoreSettings) (config.Settings, error) {
	if in == nil {
		return config.Settings{}, errors.New("core_settings must not be nil")
	}
	settings := config.Settings{
		RunID:          in.GetRunId(),
		RunDir:         in.GetRunDir(),
		SectionRule:    in.GetSectionRule(),
		RecordBatch:    in.GetRecordBatch(),
		RecordInterval: in.GetRecordInterval(),
		SaveSize:       in.GetSaveSize(),
		SaveSplit:      in.GetSaveSplit(),
		SavePart:       in.GetSavePart(),
		SaveBatch:      in.GetSaveBatch(),
		SkipStore:      in.GetSkipStore(),
		APIKey:         in.GetApiKey(),
		APIHost:        in.GetApiHost(),
		WebHost:        in.GetWebHost(),
	}
	if settings.APIHost == "" {
		return config.Settings{}, errors.New("api_host must not be empty")
	}
	if settings.APIKey == "" {
		return config.Settings{}, errors.New("api_key must not be empty")
	}
	if settings.SaveSize <= 0 || settings.SaveSplit <= 0 || settings.SaveBatch <= 0 {
		return config.Settings{}, fmt.Errorf(
			"save_size (%d), save_split (%d) and save_batch (%d) must be positive",
			settings.SaveSize, settings.SaveSplit, settings.SaveBatch)
	}
	if err := settings.Validate(); err != nil {
		return config.Settings{}, err
	}
	return settings, nil
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

// DeliverRunFinish 标记指定 session_id 的会话结束。
func (s *CoreService) DeliverRunFinish(
	_ context.Context,
	req *corev1.DeliverRunFinishRequest,
) (*corev1.DeliverRunFinishResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if err := s.runs.Finish(req.GetSessionId()); err != nil {
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
	runState, err := s.runs.State(req.GetSessionId())
	if err != nil {
		return nil, managerStatus(err)
	}
	return &corev1.GetOperationStatsResponse{
		Success: true,
		Stats:   &operationv1.OperationStats{State: coreState(runState)},
	}, nil
}

// coreState 把内部生命周期状态映射为对外的 CoreState。
// DRAINING 仍对外表现 RUNNING，直到 uploader 真正排空。
func coreState(state manager.RunState) operationv1.CoreState {
	switch state {
	case manager.RunStateDrained, manager.RunStateFailed:
		return operationv1.CoreState_CORE_STATE_FINISHED
	default:
		return operationv1.CoreState_CORE_STATE_RUNNING
	}
}

// ConfirmRunFinish 确认会话已完成并释放相关资源。
func (s *CoreService) ConfirmRunFinish(
	ctx context.Context,
	req *corev1.ConfirmRunFinishRequest,
) (*corev1.ConfirmRunFinishResponse, error) {
	if err := s.requireReady(); err != nil {
		return nil, err
	}
	if err := s.runs.Confirm(ctx, req.GetSessionId()); err != nil {
		return nil, managerStatus(err)
	}
	return &corev1.ConfirmRunFinishResponse{Success: true}, nil
}

// managerStatus 将 manager 层的领域错误转换为对应的 gRPC 状态码。
func managerStatus(err error) error {
	switch {
	case errors.Is(err, manager.ErrEmptySessionID):
		return status.Error(codes.InvalidArgument, err.Error())
	case errors.Is(err, manager.ErrSessionNotFound):
		return status.Error(codes.NotFound, err.Error())
	case errors.Is(err, manager.ErrRunNotFinished):
		return status.Error(codes.FailedPrecondition, err.Error())
	case errors.Is(err, manager.ErrRunFailed):
		return status.Error(codes.FailedPrecondition, err.Error())
	case errors.Is(err, manager.ErrInvalidSettings), errors.Is(err, manager.ErrInvalidStartRecord):
		return status.Error(codes.InvalidArgument, err.Error())
	case errors.Is(err, client.ErrMissingAPIKey), errors.Is(err, client.ErrMissingAPIHost):
		return status.Error(codes.InvalidArgument, err.Error())
	case errors.Is(err, manager.ErrTooManySessions):
		return status.Error(codes.ResourceExhausted, err.Error())
	case client.IsAuthenticationError(err):
		return status.Error(codes.Unauthenticated, err.Error())
	case client.IsNotFound(err):
		return status.Error(codes.NotFound, err.Error())
	case client.IsRetryable(err):
		return status.Error(codes.Unavailable, err.Error())
	case errors.Is(err, context.DeadlineExceeded):
		return status.Error(codes.DeadlineExceeded, err.Error())
	case errors.Is(err, context.Canceled):
		return status.Error(codes.Canceled, err.Error())
	default:
		return status.Error(codes.Internal, err.Error())
	}
}
