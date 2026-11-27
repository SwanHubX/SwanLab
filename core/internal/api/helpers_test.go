package api

import (
	"context"
	"encoding/json"
	"net/http"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
)

// verifyRoute 是 ApiKey 探活端点（api_host + /api + 业务路径）。
const verifyRoute = "/api/auth/verify"

// apiKeyProfile 是 ApiKey 模式下 /api/auth/verify 的合法响应。
func apiKeyProfile(username string) map[string]any {
	return map[string]any{
		"uid":       42,
		"username":  username,
		"createdAt": "2026-01-01T00:00:00.000Z",
	}
}

// decodeJSONBody 把记录的请求体解析成 map，便于断言字段取舍（omitempty）是否正确。
func decodeJSONBody(t *testing.T, raw []byte) map[string]any {
	t.Helper()
	var body map[string]any
	if err := json.Unmarshal(raw, &body); err != nil {
		t.Fatalf("decode request body %s: %v", raw, err)
	}
	return body
}

// newTestClient 构造一个 ApiKey 模式、带假时钟的 api 客户端。
func newTestClient(
	t *testing.T,
	cloud *clienttest.Cloud,
	opts client.Options,
) (*Client, *clienttest.Clock) {
	t.Helper()
	clock := clienttest.NewClock()
	opts.Now = clock.Now
	opts.Sleep = clock.Sleep
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	transport, err := client.New(context.Background(), cloud.URL(), "smoke-key", opts)
	if err != nil {
		t.Fatalf("create client: %v", err)
	}
	return New(transport), clock
}
