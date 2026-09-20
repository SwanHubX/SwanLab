package server

import (
	"context"
	"crypto/hmac"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/metadata"
	"google.golang.org/grpc/status"
)

// AuthTokenMetadataKey 是每个 RPC 必须携带的 auth token metadata 键，
// token 值来自 port-file 的 auth 字段。
const AuthTokenMetadataKey = "x-swanlab-core-auth-token"

// UnaryAuthInterceptor 校验每个 unary RPC 的 auth token，使用常量时间比较，
// 拒绝缺失或错误的 token。expectedToken 为空表示未启用鉴权，仅限 --listen
// 手动调试场景；SDK 启动约定（--port-filename）总会生成并配置 token。
func UnaryAuthInterceptor(expectedToken string) grpc.UnaryServerInterceptor {
	return func(
		ctx context.Context,
		req any,
		_ *grpc.UnaryServerInfo,
		handler grpc.UnaryHandler,
	) (any, error) {
		if expectedToken != "" && !requestAuthorized(ctx, expectedToken) {
			return nil, status.Error(codes.PermissionDenied, "invalid auth token")
		}
		return handler(ctx, req)
	}
}

// requestAuthorized 从 incoming metadata 中提取 auth token 并做常量时间比较。
// token 值不得出现在日志或错误消息中。
func requestAuthorized(ctx context.Context, expectedToken string) bool {
	md, ok := metadata.FromIncomingContext(ctx)
	if !ok {
		return false
	}
	values := md.Get(AuthTokenMetadataKey)
	if len(values) != 1 {
		return false
	}
	return hmac.Equal([]byte(values[0]), []byte(expectedToken))
}
