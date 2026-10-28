// Package client 是 SwanLab 云端的 HTTP 传输层：连接池复用、身份注入与底层重试。
//
// 本包只提供「一次 HTTP 往返」（Do），不含业务端点；REST 协议封装与请求/响应体
// 模型在 internal/api。
package client

import (
	"bytes"
	"context"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"io"
	"net/http"
	"net/url"
	"strings"
	"sync"
	"time"
)

const (
	// apiPrefix 是云端 REST 的固定前缀：api_host + /api + 业务路径。
	apiPrefix = "/api"
	// sdkVersionHeader 与 Python 保持一致，便于服务端按 SDK 版本分流。
	sdkVersionHeader = "X-SwanLab-SDK-Version"
	// maxResponseBytes 是单次响应体的读取上限，防止异常响应把内存打满。
	maxResponseBytes = 8 << 20
	// DefaultSDKVersion 是未经构建注入时的版本占位。
	DefaultSDKVersion = "dev"

	// 共享连接池的空闲连接上限。上传管线同时有多条指标轨与较高的文件并发，
	// http.DefaultTransport 的 2 会在持续并发下反复建连。
	defaultMaxIdleConns        = 64
	defaultMaxIdleConnsPerHost = 16
)

// Options 配置连接池、重试与可测试性注入点。零值即生产默认。
type Options struct {
	// Timeout 是单次 HTTP 请求的超时上限（含底层重试之间的等待）。
	Timeout time.Duration
	// MaxRetries 是底层 HTTP 重试次数（不含首次请求），0 表示使用 DefaultMaxRetries。
	// 需要针对单次请求关闭重试时用 NoRetries() 填入 Request.Retries，而不是把这里设成 0。
	MaxRetries int
	// BackoffFactor 是指数退避基数：第 n 次重试等待 factor * 2^(n-1)。
	BackoffFactor time.Duration
	// MaxBackoff 是单次退避等待上限。
	MaxBackoff time.Duration
	// RetryStatusCodes 覆盖默认自动重试的状态码集合。
	RetryStatusCodes []int
	// Version 是 X-SwanLab-SDK-Version 的取值，空值时使用 DefaultSDKVersion。
	Version string
	// BaseTransport 替换底层连接池（默认克隆 http.DefaultTransport）。
	BaseTransport http.RoundTripper
	// Now 注入时钟，测试用假时钟断言退避时间线。
	Now func() time.Time
	// Sleep 注入等待实现，测试可直接返回以跳过真实等待。
	Sleep func(ctx context.Context, d time.Duration) error
	// OnRetry 观测每次重试（日志或断言）。
	OnRetry func(RetryEvent)
}

// withDefaults 填充零值字段，保证 Options 可以直接零值使用。
func (o Options) withDefaults() Options {
	if o.Timeout <= 0 {
		o.Timeout = DefaultTimeout
	}
	if o.MaxRetries <= 0 {
		o.MaxRetries = DefaultMaxRetries
	}
	if o.BackoffFactor <= 0 {
		o.BackoffFactor = DefaultBackoffFactor
	}
	if o.MaxBackoff <= 0 {
		o.MaxBackoff = DefaultMaxBackoff
	}
	if len(o.RetryStatusCodes) == 0 {
		o.RetryStatusCodes = DefaultRetryStatusCodes
	}
	if o.Version == "" {
		o.Version = DefaultSDKVersion
	}
	if o.Now == nil {
		o.Now = time.Now
	}
	return o
}

// retryPolicy 把 Options 转换成底层重试策略。
func (o Options) retryPolicy() RetryPolicy {
	return RetryPolicy{
		MaxRetries:    o.MaxRetries,
		BackoffFactor: o.BackoffFactor,
		MaxBackoff:    o.MaxBackoff,
		StatusCodes:   statusSet(o.RetryStatusCodes),
		Sleep:         o.Sleep,
		OnRetry:       o.OnRetry,
	}
}

// AuthMode 是客户端实际使用的鉴权形态。
type AuthMode string

const (
	// AuthModeAPIKey 表示目标形态：常驻 Authorization: ApiKey 请求头。
	AuthModeAPIKey AuthMode = "apikey"
	// AuthModeSessionCookie 表示过渡期兜底：sid Cookie。
	AuthModeSessionCookie AuthMode = "session-cookie"
)

// Client 是一个已认证的云端客户端，绑定单一 (api_host, api_key) 身份。
//
// 并发安全：Client 自身无可变状态（身份在构造时确定），可被同账号的多个会话共享；
// 过渡期 Cookie 模式下的 sid 由其 authenticator 内部加锁维护。
type Client struct {
	apiHost string
	apiBase string
	version string
	http    *http.Client
	auth    authenticator
	now     func() time.Time
}

// New 构造并完成一次鉴权探活，返回可直接使用的 Client。
//
// 探活失败（凭证无效、后端不可达）时返回错误，不返回半初始化的客户端，
// 避免把「看起来能用」的实例交给会话层。
func New(ctx context.Context, apiHost, apiKey string, opts Options) (*Client, error) {
	opts = opts.withDefaults()
	normalized, err := normalizeHost(apiHost)
	if err != nil {
		return nil, err
	}
	if strings.TrimSpace(apiKey) == "" {
		return nil, ErrMissingAPIKey
	}
	transport := newRetryTransport(opts.BaseTransport, opts.retryPolicy())
	httpClient := &http.Client{Transport: transport, Timeout: opts.Timeout}
	auth, err := authenticate(ctx, httpClient, normalized+apiPrefix, apiKey, opts.Version, opts.Now)
	if err != nil {
		return nil, err
	}
	return &Client{
		apiHost: normalized,
		apiBase: normalized + apiPrefix,
		version: opts.Version,
		http:    httpClient,
		auth:    auth,
		now:     opts.Now,
	}, nil
}

// APIHost 返回规范化后的 api_host。
func (c *Client) APIHost() string { return c.apiHost }

// Username 返回鉴权得到的用户名。
func (c *Client) Username() string {
	if profile := c.auth.profile(); profile != nil {
		return profile.Username
	}
	return ""
}

// AuthMode 返回实际协商到的鉴权形态，供会话层记录（旧部署需提示升级）。
func (c *Client) AuthMode() AuthMode {
	if _, ok := c.auth.(*apiKeyAuth); ok {
		return AuthModeAPIKey
	}
	return AuthModeSessionCookie
}

// Now 返回客户端的时钟。
//
// 协议要求 createdAt / finishedAt 取「发起请求的时刻」，由上层（internal/api）
// 通过本方法取当前时间，保证与重试退避共用同一个可注入时钟。
func (c *Client) Now() time.Time { return c.now() }

// Request 描述一次云端调用。
//
// Path 以 "/" 开头，拼接在 api_host + /api 之后；Query 与 Body 按需填。
type Request struct {
	Method  string
	Path    string
	Query   map[string]string
	Body    any
	Retries *int // 请求级重试次数覆盖，上传类接口传 NoRetries()
}

// NoRetries 返回「关闭底层重试」的请求级覆盖（对应 Python 的 retries=0）。
//
// 上传类接口统一由上层管线的指数退避控制重试，底层不再叠加，避免请求量放大。
func NoRetries() *int {
	zero := 0
	return &zero
}

// Response 是读完并关闭响应体后的响应快照。
type Response struct {
	Method   string
	Endpoint string
	Status   int
	Header   http.Header
	Body     []byte
}

// Decode 把响应体解码到 target；空响应体视为无内容并保持 target 零值。
func (r *Response) Decode(target any) error {
	if len(bytes.TrimSpace(r.Body)) == 0 {
		return nil
	}
	if err := json.Unmarshal(r.Body, target); err != nil {
		return NewProtocolError(r.Method, r.Endpoint, "invalid_response_payload", err.Error(), err)
	}
	return nil
}

// Do 执行一次云端调用：注入身份与版本头，非 2xx 转成 APIError。
// 这是传输层唯一的请求原语，业务端点与响应语义由 internal/api 封装。
func (c *Client) Do(ctx context.Context, spec Request) (*Response, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	payload, err := encodeBody(spec.Body)
	if err != nil {
		return nil, NewProtocolError(spec.Method, spec.Path, "invalid_request_payload", err.Error(), err)
	}
	endpoint := c.apiBase + spec.Path
	if len(spec.Query) > 0 {
		endpoint += "?" + encodeQuery(spec.Query)
	}
	goCtx := ctx
	if spec.Retries != nil {
		goCtx = WithRetries(ctx, *spec.Retries)
	}

	res, err := c.send(goCtx, spec.Method, endpoint, payload)
	if err != nil {
		return nil, err
	}
	if res.Status < http.StatusOK || res.Status >= http.StatusMultipleChoices {
		return nil, newAPIError(res.Method, res.Endpoint, res.Status, res.Header, res.Body)
	}
	return res, nil
}

// send 发出单次请求并读完全部响应体。
func (c *Client) send(ctx context.Context, method, endpoint string, payload []byte) (*Response, error) {
	var body io.Reader
	if payload != nil {
		body = bytes.NewReader(payload)
	}
	req, err := http.NewRequestWithContext(ctx, method, endpoint, body)
	if err != nil {
		return nil, newTransportError(method, endpoint, err)
	}
	req.Header.Set(sdkVersionHeader, c.version)
	req.Header.Set("Accept", "application/json")
	if payload != nil {
		req.Header.Set("Content-Type", "application/json")
	}
	if authErr := c.auth.authorize(ctx, req); authErr != nil {
		return nil, authErr
	}

	resp, err := c.http.Do(req)
	if err != nil {
		// 调用方主动取消/超时优先返回 ctx 错误，保持与 Go 生态一致的中断语义。
		if ctxErr := ctx.Err(); ctxErr != nil {
			return nil, ctxErr
		}
		return nil, newTransportError(method, endpoint, err)
	}
	defer func() { _ = resp.Body.Close() }()

	raw, readErr := io.ReadAll(io.LimitReader(resp.Body, maxResponseBytes+1))
	if readErr != nil {
		if ctxErr := ctx.Err(); ctxErr != nil {
			return nil, ctxErr
		}
		return nil, newTransportError(method, endpoint, readErr)
	}
	if len(raw) > maxResponseBytes {
		return nil, NewProtocolError(method, endpoint, "response_too_large", "response body exceeds 8 MiB", nil)
	}
	return &Response{
		Method:   method,
		Endpoint: endpoint,
		Status:   resp.StatusCode,
		Header:   resp.Header,
		Body:     raw,
	}, nil
}

// encodeBody 把请求体编码为 JSON；nil 表示无请求体。
func encodeBody(body any) ([]byte, error) {
	if body == nil {
		return nil, nil
	}
	if raw, ok := body.([]byte); ok {
		return raw, nil
	}
	return json.Marshal(body)
}

// encodeQuery 按稳定顺序编码查询参数。
func encodeQuery(query map[string]string) string {
	values := url.Values{}
	for key, value := range query {
		values.Set(key, value)
	}
	return values.Encode()
}

// Factory 按 (api_host, sha256(api_key)) 复用已认证 Client，并让所有 Client
// 共享同一个底层连接池。
//
// 同账号多会话复用同一 Client 与连接；不同账号或不同站点严格隔离，杜绝串号。
type Factory struct {
	opts      Options
	transport http.RoundTripper

	mu      sync.Mutex
	clients map[string]*Client
}

// NewFactory 创建共享连接池的 Client 工厂。opts 中的 BaseTransport 决定底层
// 连接池实现，未指定时克隆 http.DefaultTransport（避免污染全局默认值）。
func NewFactory(opts Options) *Factory {
	opts = opts.withDefaults()
	base := opts.BaseTransport
	if base == nil {
		base = cloneDefaultTransport()
	}
	// 重试挂在共享连接池之上：所有 Client 共享同一份重试策略与连接复用。
	return &Factory{
		opts:      opts,
		transport: newRetryTransport(base, opts.retryPolicy()),
		clients:   make(map[string]*Client),
	}
}

// Client 返回指定身份的客户端，已缓存则直接复用。
//
// 鉴权探活在锁外完成，慢探活不会阻塞其它身份；同 key 并发可能重复探活，
// 最终只保留首个写入缓存的实例（transport 共享，丢弃不额外占用连接池）。
func (f *Factory) Client(ctx context.Context, apiHost, apiKey string) (*Client, error) {
	normalized, err := normalizeHost(apiHost)
	if err != nil {
		return nil, err
	}
	if strings.TrimSpace(apiKey) == "" {
		return nil, ErrMissingAPIKey
	}
	key := identityKey(normalized, apiKey)

	f.mu.Lock()
	existing, cached := f.clients[key]
	f.mu.Unlock()
	if cached {
		return existing, nil
	}

	httpClient := &http.Client{Transport: f.transport, Timeout: f.opts.Timeout}
	auth, err := authenticate(ctx, httpClient, normalized+apiPrefix, apiKey, f.opts.Version, f.opts.Now)
	if err != nil {
		return nil, err
	}
	created := &Client{
		apiHost: normalized,
		apiBase: normalized + apiPrefix,
		version: f.opts.Version,
		http:    httpClient,
		auth:    auth,
		now:     f.opts.Now,
	}

	f.mu.Lock()
	defer f.mu.Unlock()
	if existing, ok := f.clients[key]; ok {
		return existing, nil
	}
	f.clients[key] = created
	return created, nil
}

// Close 释放空闲连接；已取出的 Client 不因此失效，但不应在关闭后继续发起新请求。
func (f *Factory) Close() error {
	if closer, ok := f.transport.(interface{ CloseIdleConnections() }); ok {
		closer.CloseIdleConnections()
	}
	return nil
}

// identityKey 生成「站点 + 凭证摘要」的身份键：api_key 不回显在键里，避免泄漏到日志或堆栈。
func identityKey(apiHost, apiKey string) string {
	sum := sha256.Sum256([]byte(apiKey))
	return apiHost + "|" + hex.EncodeToString(sum[:])
}

// cloneDefaultTransport 克隆默认连接池，避免影响其它使用 http.DefaultTransport 的代码。
// 克隆后调高空闲连接上限，适配上传管线的并发轨道与文件上传。
func cloneDefaultTransport() *http.Transport {
	transport := &http.Transport{}
	if base, ok := http.DefaultTransport.(*http.Transport); ok {
		transport = base.Clone()
	}
	transport.MaxIdleConns = defaultMaxIdleConns
	transport.MaxIdleConnsPerHost = defaultMaxIdleConnsPerHost
	return transport
}
