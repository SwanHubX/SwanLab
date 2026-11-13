package client

import (
	"context"
	"net/http"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
)

// 测试用到的完整路由（api_host + /api + 业务路径）。
const (
	verifyRoute = apiPrefix + verifyPath
	loginRoute  = apiPrefix + loginPath
)

// apiKeyProfile 是 ApiKey 模式下 /api/auth/verify 的合法响应。
func apiKeyProfile(username string) map[string]any {
	return map[string]any{
		"uid":       42,
		"username":  username,
		"createdAt": "2026-01-01T00:00:00.000Z",
	}
}

// newAuthedClient 构造一个 ApiKey 模式、带假时钟的客户端。
func newAuthedClient(t *testing.T, cloud *clienttest.Cloud, opts Options) (*Client, *clienttest.Clock) {
	t.Helper()
	clock := clienttest.NewClock()
	opts.Now = clock.Now
	opts.Sleep = clock.Sleep
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	authed, err := New(context.Background(), cloud.URL(), "smoke-key", opts)
	if err != nil {
		t.Fatalf("create client: %v", err)
	}
	return authed, clock
}
