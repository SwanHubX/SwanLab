package client

import (
	"context"
	"errors"
	"fmt"
	"io"
	"net/http"
	"net/http/httptest"
	"strings"
	"sync"
	"sync/atomic"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
)

// metricsRequest 是传输层测试使用的通用请求：只关心「发一次 POST」，
// 业务端点语义由 internal/api 的测试覆盖。
func metricsRequest() Request {
	return Request{
		Method: http.MethodPost,
		Path:   "/house/metrics",
		Body:   map[string]any{"type": "scalar", "projectId": "proj"},
	}
}

// TestNewWithApiKey 覆盖目标形态：构造期单次 verify 探活，且后续每个请求都带
// Authorization: ApiKey 与 SDK 版本头。
func TestNewWithApiKey(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	cloud.ServeJSON(http.MethodPost, "/api/house/metrics", http.StatusOK, map[string]any{"ok": true})

	authed, clock := newAuthedClient(t, cloud, Options{})
	if got := authed.Username(); got != "alice" {
		t.Fatalf("username = %q, want alice", got)
	}
	if authed.APIHost() != cloud.URL() {
		t.Fatalf("api host = %q, want %q", authed.APIHost(), cloud.URL())
	}

	verifyReqs := cloud.RequestsFor(http.MethodGet, verifyRoute)
	if len(verifyReqs) != 1 {
		t.Fatalf("verify requests = %d, want 1 (probe must not retry)", len(verifyReqs))
	}
	if got := verifyReqs[0].Header.Get("Authorization"); got != "ApiKey smoke-key" {
		t.Fatalf("verify Authorization = %q, want %q", got, "ApiKey smoke-key")
	}
	if got := verifyReqs[0].Header.Get(sdkVersionHeader); got != DefaultSDKVersion {
		t.Fatalf("%s = %q, want %q", sdkVersionHeader, got, DefaultSDKVersion)
	}

	if _, err := authed.Do(context.Background(), metricsRequest()); err != nil {
		t.Fatalf("do request: %v", err)
	}
	reqs := cloud.RequestsFor(http.MethodPost, "/api/house/metrics")
	if len(reqs) != 1 {
		t.Fatalf("metrics requests = %d, want 1", len(reqs))
	}
	if got := reqs[0].Header.Get("Authorization"); got != "ApiKey smoke-key" {
		t.Fatalf("metrics Authorization = %q, want %q", got, "ApiKey smoke-key")
	}
	if got := reqs[0].Header.Get("Content-Type"); got != "application/json" {
		t.Fatalf("Content-Type = %q, want application/json", got)
	}
	if !strings.Contains(string(reqs[0].Body), `"projectId":"proj"`) {
		t.Fatalf("metrics body = %s, want projectId", reqs[0].Body)
	}
	if len(clock.Sleeps()) != 0 {
		t.Fatalf("unexpected retries: %v", clock.Sleeps())
	}
}

// TestNewRejectsInvalidCredentials 覆盖 401：凭证被拒绝必须阻断构造，
// 不返回半初始化的客户端。
func TestNewRejectsInvalidCredentials(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusUnauthorized, map[string]string{
		"code": "unauthorized", "message": "invalid api key",
	})

	authed, err := New(context.Background(), cloud.URL(), "bad-key", Options{Now: clienttest.NewClock().Now})
	if authed != nil {
		t.Fatal("client must be nil when authentication fails")
	}
	if !IsAuthenticationError(err) {
		t.Fatalf("error = %v, want authentication error", err)
	}
	if IsRetryable(err) {
		t.Fatal("authentication failure must not be retryable")
	}
	if !strings.Contains(err.Error(), "invalid api key") {
		t.Fatalf("error = %v, want backend message preserved", err)
	}
}

// TestNewRejectsInvalidProfile 覆盖「响应结构不合法」：uid 非正数属于协议错误，
// 不能建立运行时客户端。
func TestNewRejectsInvalidProfile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, map[string]any{
		"uid": 0, "username": "alice", "createdAt": "2026-01-01T00:00:00.000Z",
	})

	_, err := New(context.Background(), cloud.URL(), "smoke-key", Options{})
	if !IsInvalidResponse(err) {
		t.Fatalf("error = %v, want invalid response", err)
	}
	if IsRetryable(err) {
		t.Fatal("protocol error must not be retryable")
	}
}

// TestLegacyCookieFallback 覆盖过渡期兜底：后端缺少 /api/auth/verify 时退化为
// sid Cookie 登录，后续请求携带 Cookie 而不是 ApiKey 头。
func TestLegacyCookieFallback(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, loginRoute, http.StatusOK, map[string]any{
		"sid":       "sid-xyz",
		"expiredAt": "2999-01-01T00:00:00.000000Z",
		"userInfo":  map[string]string{"username": "legacy"},
	})
	cloud.ServeJSON(http.MethodPost, "/api/house/metrics", http.StatusOK, map[string]any{"ok": true})

	clock := clienttest.NewClock()
	authed, err := New(context.Background(), cloud.URL(), "legacy-key", Options{Now: clock.Now, Sleep: clock.Sleep})
	if err != nil {
		t.Fatalf("create legacy client: %v", err)
	}
	if got := authed.Username(); got != "legacy" {
		t.Fatalf("username = %q, want legacy", got)
	}
	if got := authed.AuthMode(); got != AuthModeSessionCookie {
		t.Fatalf("auth mode = %q, want %q", got, AuthModeSessionCookie)
	}

	loginReqs := cloud.RequestsFor(http.MethodPost, loginRoute)
	if len(loginReqs) != 1 {
		t.Fatalf("login requests = %d, want 1", len(loginReqs))
	}
	if got := loginReqs[0].Header.Get("authorization"); got != "legacy-key" {
		t.Fatalf("login authorization = %q, want raw api key", got)
	}

	if _, err := authed.Do(context.Background(), metricsRequest()); err != nil {
		t.Fatalf("do request: %v", err)
	}
	metricsReqs := cloud.RequestsFor(http.MethodPost, "/api/house/metrics")
	if len(metricsReqs) != 1 {
		t.Fatalf("metrics requests = %d, want 1", len(metricsReqs))
	}
	if got := metricsReqs[0].Header.Get("Authorization"); got != "" {
		t.Fatalf("legacy mode must not send ApiKey header, got %q", got)
	}
	if cookie := metricsReqs[0].Header.Get("Cookie"); !strings.Contains(cookie, "sid=sid-xyz") {
		t.Fatalf("legacy Cookie = %q, want sid cookie", cookie)
	}
}

// TestLegacyCookieDoesNotReloginOnUnauthorized 覆盖过渡期 401：与 Python 一致，
// 请求层不重放，直接返回鉴权错误，也不会额外登录。
func TestLegacyCookieDoesNotReloginOnUnauthorized(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var logins atomic.Int32
	cloud.Handle(http.MethodPost, loginRoute, func(w http.ResponseWriter, _ *http.Request) {
		n := logins.Add(1)
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{
			"sid":       fmt.Sprintf("sid-%d", n),
			"expiredAt": "2999-01-01T00:00:00.000000Z",
			"userInfo":  map[string]string{"username": "legacy"},
		})
	})
	var hits atomic.Int32
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, _ *http.Request) {
		hits.Add(1)
		clienttest.WriteJSON(w, http.StatusUnauthorized, map[string]string{"code": "unauthorized", "message": "sid expired"})
	})

	clock := clienttest.NewClock()
	authed, err := New(context.Background(), cloud.URL(), "legacy-key", Options{Now: clock.Now, Sleep: clock.Sleep})
	if err != nil {
		t.Fatalf("create legacy client: %v", err)
	}
	_, err = authed.Do(context.Background(), metricsRequest())
	if !IsAuthenticationError(err) {
		t.Fatalf("error = %v, want authentication error", err)
	}
	if got := hits.Load(); got != 1 {
		t.Fatalf("metrics attempts = %d, want 1 (401 must not be replayed)", got)
	}
	if got := logins.Load(); got != 1 {
		t.Fatalf("logins = %d, want 1 (initial login only)", got)
	}
}

// TestMissingCredentials 覆盖入参校验。
func TestMissingCredentials(t *testing.T) {
	if _, err := New(context.Background(), "https://api.example.com", "", Options{}); !errors.Is(err, ErrMissingAPIKey) {
		t.Fatalf("error = %v, want ErrMissingAPIKey", err)
	}
	if _, err := New(context.Background(), "  ", "key", Options{}); !errors.Is(err, ErrMissingAPIHost) {
		t.Fatalf("error = %v, want ErrMissingAPIHost", err)
	}
}

// TestRetryOnServiceUnavailable 覆盖「503 → 200」：底层自动重试并遵循 0.5 起步的
// 指数退避，最终请求成功。
func TestRetryOnServiceUnavailable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	var attempts atomic.Int32
	cloud.Handle(http.MethodPost, "/api/house/experiments/exp-1/heartbeat", func(w http.ResponseWriter, _ *http.Request) {
		if attempts.Add(1) <= 2 {
			clienttest.WriteJSON(w, http.StatusServiceUnavailable, map[string]string{"code": "unavailable"})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	var events []RetryEvent
	authed, clock := newAuthedClient(t, cloud, Options{
		OnRetry: func(e RetryEvent) { events = append(events, e) },
	})
	_, err := authed.Do(context.Background(), Request{
		Method: http.MethodPost,
		Path:   "/house/experiments/exp-1/heartbeat",
		Body:   map[string]string{"flagId": "123456"},
	})
	if err != nil {
		t.Fatalf("do request: %v", err)
	}
	if got := cloud.CountFor(http.MethodPost, "/api/house/experiments/exp-1/heartbeat"); got != 3 {
		t.Fatalf("attempts = %d, want 3", got)
	}
	want := []time.Duration{500 * time.Millisecond, time.Second}
	if got := clock.Sleeps(); len(got) != len(want) || got[0] != want[0] || got[1] != want[1] {
		t.Fatalf("backoff = %v, want %v", got, want)
	}
	if len(events) != 2 || events[0].Status != http.StatusServiceUnavailable || events[0].Attempt != 1 {
		t.Fatalf("retry events = %+v, want two 503 events", events)
	}
}

// TestRetryExhaustedKeepsRetryable 覆盖「重试耗尽」：底层停止重试后，管线层仍应看到
// 可重试错误，从而保留记录。
func TestRetryExhaustedKeepsRetryable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	cloud.ServeJSON(http.MethodPost, "/api/house/experiments/exp-1/heartbeat", http.StatusServiceUnavailable, map[string]string{"code": "unavailable"})

	authed, clock := newAuthedClient(t, cloud, Options{MaxRetries: 2})
	_, err := authed.Do(context.Background(), Request{
		Method: http.MethodPost,
		Path:   "/house/experiments/exp-1/heartbeat",
	})
	if err == nil {
		t.Fatal("expected error")
	}
	if !IsRetryable(err) {
		t.Fatalf("error = %v, want retryable", err)
	}
	if got := cloud.CountFor(http.MethodPost, "/api/house/experiments/exp-1/heartbeat"); got != 3 {
		t.Fatalf("attempts = %d, want 3 (1 + 2 retries)", got)
	}
	if got := len(clock.Sleeps()); got != 2 {
		t.Fatalf("sleeps = %d, want 2", got)
	}
}

// TestNoRetriesOverrideKeepsRetryable 覆盖「请求级关闭底层重试」：上传类接口把重试权
// 交给管线层，避免双层重试叠加放大请求量（对齐 Python retries=0）。
func TestNoRetriesOverrideKeepsRetryable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	cloud.ServeJSON(http.MethodPost, "/api/house/metrics", http.StatusServiceUnavailable, map[string]string{"code": "unavailable"})

	authed, clock := newAuthedClient(t, cloud, Options{})
	spec := metricsRequest()
	spec.Retries = NoRetries()
	_, err := authed.Do(context.Background(), spec)
	if !IsRetryable(err) {
		t.Fatalf("error = %v, want retryable", err)
	}
	if got := cloud.CountFor(http.MethodPost, "/api/house/metrics"); got != 1 {
		t.Fatalf("attempts = %d, want 1 (no transport retry)", got)
	}
	if got := len(clock.Sleeps()); got != 0 {
		t.Fatalf("sleeps = %v, want none", clock.Sleeps())
	}
}

// TestRetryAfterHonored 覆盖服务端冷却窗口：Retry-After 覆盖指数退避。
func TestRetryAfterHonored(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	var attempts atomic.Int32
	cloud.Handle(http.MethodPost, "/api/house/experiments/exp-1/heartbeat", func(w http.ResponseWriter, _ *http.Request) {
		if attempts.Add(1) == 1 {
			w.Header().Set("Retry-After", "3")
			clienttest.WriteJSON(w, http.StatusTooManyRequests, map[string]string{"code": "rate_limited"})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	authed, clock := newAuthedClient(t, cloud, Options{})
	if _, err := authed.Do(context.Background(), Request{
		Method: http.MethodPost,
		Path:   "/house/experiments/exp-1/heartbeat",
	}); err != nil {
		t.Fatalf("do request: %v", err)
	}
	if got := clock.Sleeps(); len(got) != 1 || got[0] != 3*time.Second {
		t.Fatalf("sleeps = %v, want [3s] from Retry-After", got)
	}
}

// TestFactoryReusesAndIsolates 覆盖身份隔离：同 (host, key) 复用同一 Client，
// 不同 key 不共享 Client，但共享同一个底层连接池。
func TestFactoryReusesAndIsolates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))

	factory := NewFactory(Options{Now: clienttest.NewClock().Now})
	t.Cleanup(func() { _ = factory.Close() })

	first, err := factory.Client(context.Background(), cloud.URL(), "key-a")
	if err != nil {
		t.Fatalf("acquire first client: %v", err)
	}
	again, err := factory.Client(context.Background(), cloud.URL(), "key-a")
	if err != nil {
		t.Fatalf("acquire reused client: %v", err)
	}
	if first != again {
		t.Fatal("same identity must reuse the same client")
	}
	other, err := factory.Client(context.Background(), cloud.URL(), "key-b")
	if err != nil {
		t.Fatalf("acquire other client: %v", err)
	}
	if other == first {
		t.Fatal("different api key must not reuse the client")
	}
	if other.http.Transport != first.http.Transport {
		t.Fatal("clients must share the process-wide transport")
	}
	// 同 host 不同 key 各自完成一次探活。
	if got := cloud.CountFor(http.MethodGet, verifyRoute); got != 2 {
		t.Fatalf("verify probes = %d, want 2 (one per identity)", got)
	}
}

// TestFactoryConcurrentSameIdentity 覆盖并发取用：同 key 的并发调用收敛到同一实例。
func TestFactoryConcurrentSameIdentity(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))

	factory := NewFactory(Options{Now: clienttest.NewClock().Now})
	t.Cleanup(func() { _ = factory.Close() })

	const callers = 8
	clients := make([]*Client, callers)
	errs := make([]error, callers)
	start := make(chan struct{})
	var wg sync.WaitGroup
	for i := 0; i < callers; i++ {
		wg.Add(1)
		go func(index int) {
			defer wg.Done()
			<-start
			clients[index], errs[index] = factory.Client(context.Background(), cloud.URL(), "key-a")
		}(i)
	}
	close(start)
	wg.Wait()

	for i := range clients {
		if errs[i] != nil {
			t.Fatalf("caller %d: %v", i, errs[i])
		}
		if clients[i] != clients[0] {
			t.Fatalf("caller %d got a different client; concurrent callers must converge on one instance", i)
		}
	}
}

// TestCloneDefaultTransportPools 校验池参数已调高，且不污染 http.DefaultTransport。
func TestCloneDefaultTransportPools(t *testing.T) {
	transport := cloneDefaultTransport()
	if got := transport.MaxIdleConnsPerHost; got != defaultMaxIdleConnsPerHost {
		t.Fatalf("MaxIdleConnsPerHost = %d, want %d", got, defaultMaxIdleConnsPerHost)
	}
	if got := transport.MaxIdleConns; got != defaultMaxIdleConns {
		t.Fatalf("MaxIdleConns = %d, want %d", got, defaultMaxIdleConns)
	}
	base, ok := http.DefaultTransport.(*http.Transport)
	if ok && base.MaxIdleConnsPerHost == defaultMaxIdleConnsPerHost {
		t.Fatal("clone must not mutate http.DefaultTransport")
	}
}

// closeTracker 是只统计 CloseIdleConnections 调用的 RoundTripper。
type closeTracker struct {
	closed atomic.Int32
}

func (c *closeTracker) RoundTrip(*http.Request) (*http.Response, error) {
	return nil, errors.New("transport is not used by this test")
}

func (c *closeTracker) CloseIdleConnections() { c.closed.Add(1) }

// TestFactoryCloseReleasesIdleConnections 覆盖 Factory.Close 会透传到共享连接池。
func TestFactoryCloseReleasesIdleConnections(t *testing.T) {
	tracker := &closeTracker{}
	factory := NewFactory(Options{BaseTransport: tracker})
	if err := factory.Close(); err != nil {
		t.Fatalf("close factory: %v", err)
	}
	if got := tracker.closed.Load(); got != 1 {
		t.Fatalf("CloseIdleConnections calls = %d, want 1", got)
	}
}

// truncatedBodyTransport 返回 2xx 但读取响应体时中断，用于覆盖鉴权链路的读错误分类。
type truncatedBodyTransport struct {
	verifyStatus int
}

func (t truncatedBodyTransport) RoundTrip(req *http.Request) (*http.Response, error) {
	status := http.StatusOK
	if strings.HasSuffix(req.URL.Path, verifyPath) {
		status = t.verifyStatus
	}
	return &http.Response{
		StatusCode: status,
		Header:     http.Header{},
		Body:       io.NopCloser(errReader{}),
		Request:    req,
	}, nil
}

type errReader struct{}

func (errReader) Read([]byte) (int, error) { return 0, io.ErrUnexpectedEOF }

// TestVerifyReadErrorIsRetryable 覆盖探活响应读取中断：按传输错误处理，保留重试。
func TestVerifyReadErrorIsRetryable(t *testing.T) {
	_, err := New(context.Background(), "https://api.example.com", "smoke-key",
		Options{BaseTransport: truncatedBodyTransport{verifyStatus: http.StatusOK}})
	if !IsRetryable(err) || IsInvalidResponse(err) {
		t.Fatalf("error = %v, want retryable transport error", err)
	}
}

// TestLoginReadErrorIsRetryable 覆盖旧后端登录响应读取中断：同样保留重试。
func TestLoginReadErrorIsRetryable(t *testing.T) {
	_, err := New(context.Background(), "https://api.example.com", "smoke-key",
		Options{BaseTransport: truncatedBodyTransport{verifyStatus: http.StatusNotFound}})
	if !IsRetryable(err) || IsInvalidResponse(err) {
		t.Fatalf("error = %v, want retryable transport error", err)
	}
}

// TestNormalizeHost 覆盖 host 归一化：补协议、剥离路由/查询串、保留端口。
func TestNormalizeHost(t *testing.T) {
	cases := []struct {
		in   string
		want string
	}{
		{"https://api.swanlab.cn", "https://api.swanlab.cn"},
		{"https://api.swanlab.cn/", "https://api.swanlab.cn"},
		{"api.swanlab.cn", "https://api.swanlab.cn"},
		{"https://api.swanlab.cn/api/v1?x=1", "https://api.swanlab.cn"},
		{"http://127.0.0.1:8080", "http://127.0.0.1:8080"},
		{"http://127.0.0.1:8080/prefix", "http://127.0.0.1:8080"},
		{"http://[::1]:8080", "http://[::1]:8080"},
		{"[::1]:8080", "https://[::1]:8080"},
		{"  https://api.swanlab.cn  ", "https://api.swanlab.cn"},
	}
	for _, tc := range cases {
		got, err := normalizeHost(tc.in)
		if err != nil {
			t.Fatalf("normalizeHost(%q): %v", tc.in, err)
		}
		if got != tc.want {
			t.Fatalf("normalizeHost(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

// TestIdentityKeyHidesSecret 覆盖身份键不含明文 api_key。
func TestIdentityKeyHidesSecret(t *testing.T) {
	key := identityKey("https://api.swanlab.cn", "super-secret-key")
	if strings.Contains(key, "super-secret-key") {
		t.Fatalf("identity key leaks api key: %q", key)
	}
	if key == identityKey("https://api.swanlab.cn", "another-key") {
		t.Fatal("different api keys must produce different identity keys")
	}
	if identityKey("https://a.example.com", "k") == identityKey("https://b.example.com", "k") {
		t.Fatal("different hosts must produce different identity keys")
	}
}

// TestContextCancellationIsNotRetried 覆盖调用方取消：不重试且返回 ctx 错误。
// 取消时机由 handler 信号驱动，避免 sleep 与请求发出竞速导致的偶发失败。
func TestContextCancellationIsNotRetried(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	requested := make(chan struct{}, 1)
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(_ http.ResponseWriter, r *http.Request) {
		select {
		case requested <- struct{}{}:
		default:
		}
		<-r.Context().Done()
	})

	authed, clock := newAuthedClient(t, cloud, Options{})
	ctx, cancel := context.WithCancel(context.Background())
	// done 在请求未发出时叫停等待协程，避免永久阻塞。
	done := make(chan struct{})
	t.Cleanup(func() {
		close(done)
		cancel()
	})
	go func() {
		select {
		case <-requested:
		case <-done:
		}
		cancel()
	}()

	_, err := authed.Do(ctx, metricsRequest())
	if !errors.Is(err, context.Canceled) {
		t.Fatalf("error = %v, want context.Canceled", err)
	}
	if got := cloud.CountFor(http.MethodPost, "/api/house/metrics"); got != 1 {
		t.Fatalf("attempts = %d, want 1", got)
	}
	if got := len(clock.Sleeps()); got != 0 {
		t.Fatalf("sleeps = %v, want none", clock.Sleeps())
	}
}

// TestApiKeyModeDoesNotReloginOnUnauthorized 覆盖 ApiKey 模式的 401：无状态凭证不重登，
// 直接返回鉴权错误且只发一次请求（过渡期 Cookie 模式才需要换发 sid）。
func TestApiKeyModeDoesNotReloginOnUnauthorized(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))
	cloud.ServeJSON(http.MethodPost, "/api/house/experiments/exp-1/heartbeat", http.StatusUnauthorized, map[string]string{
		"code": "unauthorized", "message": "api key revoked",
	})

	authed, _ := newAuthedClient(t, cloud, Options{})
	if got := authed.AuthMode(); got != AuthModeAPIKey {
		t.Fatalf("auth mode = %q, want %q", got, AuthModeAPIKey)
	}
	_, err := authed.Do(context.Background(), Request{
		Method: http.MethodPost,
		Path:   "/house/experiments/exp-1/heartbeat",
	})
	if !IsAuthenticationError(err) {
		t.Fatalf("error = %v, want authentication error", err)
	}
	if got := cloud.CountFor(http.MethodPost, "/api/house/experiments/exp-1/heartbeat"); got != 1 {
		t.Fatalf("attempts = %d, want 1 (no relogin for stateless api key)", got)
	}
	if got := cloud.CountFor(http.MethodPost, loginRoute); got != 0 {
		t.Fatalf("login attempts = %d, want 0", got)
	}
}

// TestNetworkErrorIsRetryable 覆盖连接失败：无响应时保留重试语义与 network_error 业务码，
// 供管线层退避后重新投递。
func TestNetworkErrorIsRetryable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, verifyRoute, http.StatusOK, apiKeyProfile("alice"))

	authed, _ := newAuthedClient(t, cloud, Options{})
	cloud.Close() // 关闭后端，后续请求将连接失败

	_, err := authed.Do(context.Background(), metricsRequest())
	if err == nil {
		t.Fatal("expected network error")
	}
	if StatusCode(err) != 0 {
		t.Fatalf("status = %d, want 0 for transport error", StatusCode(err))
	}
	if !IsRetryable(err) {
		t.Fatalf("error = %v, want retryable", err)
	}
	var apiErr *APIError
	if !errors.As(err, &apiErr) || apiErr.Code != "network_error" {
		t.Fatalf("error = %v, want code network_error", err)
	}
}

// TestRetryTransportClosesBodyOnRetry 确认重试会读走并关闭上一次响应体，
// 使连接可以回到池里复用（httptest 下以请求计数间接验证）。
func TestRetryTransportClosesBodyOnRetry(t *testing.T) {
	var attempts atomic.Int32
	server := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, _ *http.Request) {
		if attempts.Add(1) == 1 {
			w.WriteHeader(http.StatusBadGateway)
			_, _ = w.Write([]byte(strings.Repeat("x", 4096)))
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	}))
	t.Cleanup(server.Close)

	clock := clienttest.NewClock()
	transport := newRetryTransport(server.Client().Transport, RetryPolicy{
		MaxRetries:    1,
		BackoffFactor: 10 * time.Millisecond,
		Sleep:         clock.Sleep,
	})
	httpClient := &http.Client{Transport: transport}
	resp, err := httpClient.Get(server.URL)
	if err != nil {
		t.Fatalf("get: %v", err)
	}
	defer func() { _ = resp.Body.Close() }()
	if resp.StatusCode != http.StatusOK {
		t.Fatalf("status = %d, want 200", resp.StatusCode)
	}
	if got := attempts.Load(); got != 2 {
		t.Fatalf("attempts = %d, want 2", got)
	}
}
