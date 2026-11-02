package client

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
	"net/url"
	"strings"
	"sync"
	"time"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// 鉴权相关端点，均相对于 api_host + /api 前缀。
const (
	verifyPath = "/auth/verify"
	loginPath  = "/login/api_key"
)

// refreshBuffer 是 Cookie 模式的提前刷新窗口：距过期不足 7 天即换发 sid，
// 与 Python Client.REFRESH_TIME 保持一致。
const refreshBuffer = 7 * 24 * time.Hour

// Profile 是后端返回的当前身份信息。
//
// 只强校验鉴权必需字段（uid 为正、username 非空、createdAt 非空）；
// name / avatar / verified 允许后端省略，缺失时归一化为零值。
type Profile struct {
	UID       int64  `json:"uid"`
	Username  string `json:"username"`
	Name      string `json:"name"`
	Avatar    string `json:"avatar"`
	Verified  bool   `json:"verified"`
	CreatedAt string `json:"createdAt"`
}

// validate 校验鉴权必需字段，缺失或非法都视为协议错误（不能建立运行时客户端）。
func (p *Profile) validate(method, path string) error {
	switch {
	case p == nil:
		return NewProtocolError(method, path, "invalid_auth_profile", "empty auth profile", nil)
	case p.UID <= 0:
		return NewProtocolError(method, path, "invalid_auth_profile", "auth profile uid must be a positive integer", nil)
	case p.Username == "":
		return NewProtocolError(method, path, "invalid_auth_profile", "auth profile username must not be empty", nil)
	case p.CreatedAt == "":
		return NewProtocolError(method, path, "invalid_auth_profile", "auth profile createdAt must be a string", nil)
	}
	return nil
}

// authenticator 负责把身份注入每个请求。
//
// 两种实现的切换完全收敛在本文件：上层业务代码只认 Client，不感知鉴权形态。
// 与 Python 一致，401 不在请求层重放：凭证失效由上层按鉴权错误处理。
type authenticator interface {
	// authorize 在请求发出前注入凭证。
	authorize(ctx context.Context, req *http.Request) error
	// profile 返回初始化时校验通过的身份信息。
	profile() *Profile
}

// apiKeyAuth 是目标形态：无状态常驻请求头 Authorization: ApiKey <key>。
// ApiKey 无生命周期终结语义，因此没有刷新逻辑与 401 重登。
type apiKeyAuth struct {
	key  string
	user *Profile
}

var _ authenticator = (*apiKeyAuth)(nil)

func (a *apiKeyAuth) authorize(_ context.Context, req *http.Request) error {
	req.Header.Set("Authorization", "ApiKey "+a.key)
	return nil
}

func (a *apiKeyAuth) profile() *Profile { return a.user }

// cookieAuth 是过渡期兜底：后端尚未支持 ApiKey 请求头直连时，退化为
// POST /api/login/api_key 换取 sid Cookie，并在距过期 7 天内提前换发。
//
// 与 Python 一致，401 不触发重登；过渡完成后本实现可直接删除，上层零改动。
type cookieAuth struct {
	apiBase string
	apiKey  string
	version string
	http    *http.Client
	now     func() time.Time

	mu        sync.Mutex
	sid       string
	expiredAt time.Time
	user      *Profile
}

var _ authenticator = (*cookieAuth)(nil)

// ensureFresh 先在短临界区内判断是否临近过期，只有需要换发时才进入 login 的串行化路径。
func (a *cookieAuth) ensureFresh(ctx context.Context) error {
	a.mu.Lock()
	fresh := a.sid != "" && a.now().Add(refreshBuffer).Before(a.expiredAt)
	a.mu.Unlock()
	if fresh {
		return nil
	}
	_, err := a.login(ctx)
	return err
}

func (a *cookieAuth) authorize(ctx context.Context, req *http.Request) error {
	if err := a.ensureFresh(ctx); err != nil {
		return err
	}
	a.mu.Lock()
	sid := a.sid
	a.mu.Unlock()
	req.AddCookie(&http.Cookie{Name: "sid", Value: sid})
	return nil
}

func (a *cookieAuth) profile() *Profile {
	a.mu.Lock()
	defer a.mu.Unlock()
	return a.user
}

// login 换发 sid，串行化在 a.mu 内。
func (a *cookieAuth) login(ctx context.Context) (bool, error) {
	a.mu.Lock()
	defer a.mu.Unlock()

	endpoint := a.apiBase + loginPath
	req, err := http.NewRequestWithContext(ctx, http.MethodPost, endpoint, nil)
	if err != nil {
		return false, newTransportError(http.MethodPost, endpoint, err)
	}
	// 登录接口使用裸 api_key（与 Python bootstrap.login_by_api_key 一致）。
	req.Header.Set("authorization", a.apiKey)
	req.Header.Set(sdkVersionHeader, a.version)

	resp, err := a.http.Do(req)
	if err != nil {
		return false, newTransportError(http.MethodPost, endpoint, err)
	}
	defer func() { _ = resp.Body.Close() }()
	body, readErr := io.ReadAll(io.LimitReader(resp.Body, maxResponseBytes))
	if resp.StatusCode != http.StatusOK {
		return false, newAPIError(http.MethodPost, endpoint, resp.StatusCode, resp.Header, body)
	}
	if readErr != nil {
		return false, newTransportError(http.MethodPost, endpoint, readErr)
	}

	var payload loginResponse
	if decodeErr := json.Unmarshal(body, &payload); decodeErr != nil {
		return false, NewProtocolError(http.MethodPost, endpoint, "invalid_login_response", decodeErr.Error(), decodeErr)
	}
	if payload.SID == "" {
		return false, NewProtocolError(http.MethodPost, endpoint, "invalid_login_response", "login response missing sid", nil)
	}
	expiredAt, expireErr := parseExpiredAt(payload.ExpiredAt)
	if expireErr != nil {
		return false, NewProtocolError(http.MethodPost, endpoint, "invalid_login_response", expireErr.Error(), expireErr)
	}
	user := &Profile{Username: payload.UserInfo.Username}
	if user.Username == "" {
		return false, NewProtocolError(
			http.MethodPost, endpoint, "invalid_login_response", "login response missing userInfo.username", nil)
	}

	a.sid = payload.SID
	a.expiredAt = expiredAt
	a.user = user
	console.Debug("refreshed SwanLab session cookie")
	return true, nil
}

// loginResponse 是 POST /api/login/api_key 的响应体。
type loginResponse struct {
	SID       string `json:"sid"`
	ExpiredAt string `json:"expiredAt"`
	UserInfo  struct {
		Username string `json:"username"`
	} `json:"userInfo"`
}

// parseExpiredAt 解析后端返回的过期时间（ISO 8601，通常带微秒）。
// 兼容带/不带小数秒、带/不带 Z 的写法，避免因后端格式微调导致整体鉴权失败。
func parseExpiredAt(value string) (time.Time, error) {
	if value == "" {
		return time.Time{}, errors.New("login response missing expiredAt")
	}
	layouts := []string{
		"2006-01-02T15:04:05.999999999Z07:00",
		time.RFC3339Nano,
		"2006-01-02T15:04:05.999999999",
		"2006-01-02T15:04:05",
	}
	for _, layout := range layouts {
		if at, err := time.Parse(layout, value); err == nil {
			return at.UTC(), nil
		}
	}
	return time.Time{}, fmt.Errorf("invalid expiredAt %q", value)
}

// authenticate 探测后端鉴权能力并构造身份注入器。
//
// 判定顺序：
//  1. GET /api/auth/verify（单次请求，不重试）：200 且身份字段合法 → ApiKey 模式；
//  2. 404 说明后端太旧（缺该端点）→ 记录明确的升级提示后退化为 Cookie 模式；
//  3. 401/403 → 返回 ErrNotAuthenticated，由调用方决定是否阻断启动；
//  4. 其它错误原样返回，保留可重试语义。
func authenticate(
	ctx context.Context,
	httpClient *http.Client,
	apiBase, apiKey, version string,
	now func() time.Time,
) (authenticator, error) {
	endpoint := apiBase + verifyPath
	req, err := http.NewRequestWithContext(WithRetries(ctx, 0), http.MethodGet, endpoint, nil)
	if err != nil {
		return nil, newTransportError(http.MethodGet, endpoint, err)
	}
	req.Header.Set("Authorization", "ApiKey "+apiKey)
	req.Header.Set(sdkVersionHeader, version)

	resp, err := httpClient.Do(req)
	if err != nil {
		return nil, newTransportError(http.MethodGet, endpoint, err)
	}
	defer func() { _ = resp.Body.Close() }()
	body, readErr := io.ReadAll(io.LimitReader(resp.Body, maxResponseBytes))

	if resp.StatusCode == http.StatusNotFound {
		console.Warning(
			"backend does not expose /api/auth/verify, falling back to legacy session-cookie authentication;",
			"please upgrade your SwanLab deployment to enable ApiKey authentication.",
		)
		cookie := &cookieAuth{
			apiBase: apiBase,
			apiKey:  apiKey,
			version: version,
			http:    httpClient,
			now:     now,
		}
		if _, err := cookie.login(WithRetries(ctx, DefaultMaxRetries)); err != nil {
			return nil, err
		}
		return cookie, nil
	}
	// 读取中断属于传输层问题，保留重试语义，不能当成协议错误终结。
	if readErr != nil {
		return nil, newTransportError(http.MethodGet, endpoint, readErr)
	}
	if resp.StatusCode != http.StatusOK {
		return nil, newAPIError(http.MethodGet, endpoint, resp.StatusCode, resp.Header, body)
	}

	var profile Profile
	if decodeErr := json.Unmarshal(body, &profile); decodeErr != nil {
		return nil, NewProtocolError(http.MethodGet, endpoint, "invalid_auth_profile", decodeErr.Error(), decodeErr)
	}
	if validateErr := profile.validate(http.MethodGet, endpoint); validateErr != nil {
		return nil, validateErr
	}
	return &apiKeyAuth{key: apiKey, user: &profile}, nil
}

// normalizeHost 规范化 api_host：剥离路由与查询串，缺协议时补 https，去掉末尾斜杠。
// 与 Python nrc.fmt 保持一致，确保同一个站点不会因书写差异被当成两个身份。
func normalizeHost(raw string) (string, error) {
	host := strings.TrimSpace(raw)
	if host == "" {
		return "", ErrMissingAPIHost
	}
	candidate := host
	if !hasScheme(host) {
		candidate = "https://" + host
	}
	parsed, err := url.Parse(candidate)
	if err != nil || parsed.Host == "" {
		return "", fmt.Errorf("%w: %q", ErrMissingAPIHost, raw)
	}
	scheme := parsed.Scheme
	if scheme == "" {
		scheme = "https"
	}
	// 直接使用 parsed.Host：它保留 IPv6 的方括号与端口，且不含 userinfo。
	return scheme + "://" + parsed.Host, nil
}

// hasScheme 判断是否已经带协议前缀。
func hasScheme(host string) bool {
	return strings.Contains(host, "://")
}
