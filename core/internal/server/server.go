// Package server 提供 swanlab-core 的 gRPC 服务端实现。
//
// 当前为脚手架阶段：CoreService v2 已注册并可完成服务级关闭，
// 鉴权 interceptor、capability 应答与 run 会话路由在后续迭代接入。
package server

import (
	"context"
	"crypto/hmac"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/status"

	corev2 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v2"
)

// Service 实现 CoreService v2。
//
// 除 TeardownService 外的所有 RPC 暂由嵌入的 Unimplemented 实现接管，
// 返回 UNIMPLEMENTED；run 级接口在会话路由迭代中逐步补齐。
type Service struct {
	corev2.UnimplementedCoreServiceServer

	ownerToken string
	controller *Controller
}

// NewService 创建服务实例。ownerToken 为 spawn owner 通过私有文件传入的
// 服务级令牌，是唯一允许触发 TeardownService 的凭证。
func NewService(ownerToken string, controller *Controller) *Service {
	return &Service{
		ownerToken: ownerToken,
		controller: controller,
	}
}

// Register 将服务注册到 gRPC Server。
func (s *Service) Register(g *grpc.Server) {
	corev2.RegisterCoreServiceServer(g, s)
}

// TeardownService 关闭整个服务进程，仅接受正确的 owner token。
// 校验使用常量时间比较；token 值不得出现在日志或错误消息中。
func (s *Service) TeardownService(_ context.Context, req *corev2.TeardownServiceRequest) (*corev2.TeardownServiceResponse, error) {
	if s.ownerToken == "" || !hmac.Equal([]byte(req.OwnerToken), []byte(s.ownerToken)) {
		return nil, status.Error(codes.PermissionDenied, "invalid owner token")
	}
	// 异步触发统一关闭路径，保证本响应先于连接关闭送达调用方。
	s.controller.Shutdown("teardown")
	return &corev2.TeardownServiceResponse{}, nil
}
