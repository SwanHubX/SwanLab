package client

import (
	"context"
	"errors"
	"fmt"
	"io"
	"net/http"
	"time"
)

// 与 Python urllib3 Retry 对齐的默认参数。
const (
	// DefaultMaxRetries 是底层 HTTP 的重试次数：首次请求 + 5 次重放。
	DefaultMaxRetries = 5
	// DefaultBackoffFactor 是退避因子，第 n 次重试等待 factor * 2^(n-1)：
	// 0.5s、1s、2s、4s、8s。
	DefaultBackoffFactor = 500 * time.Millisecond
	// DefaultMaxBackoff 是单次等待上限，对齐 urllib3 的 backoff_max。
	DefaultMaxBackoff = 120 * time.Second
	// DefaultTimeout 是单次 HTTP 请求的超时上限（含底层重试之间的等待）。
	DefaultTimeout = 60 * time.Second
)

// DefaultRetryStatusCodes 是需要自动重试的服务端状态码，取自 urllib3 status_forcelist。
var DefaultRetryStatusCodes = []int{
	http.StatusTooManyRequests,
	http.StatusInternalServerError,
	http.StatusBadGateway,
	http.StatusServiceUnavailable,
	http.StatusGatewayTimeout,
}

// RetryEvent 描述一次即将发生的重试，供调用方观测（日志、测试断言）。
type RetryEvent struct {
	Method  string
	URL     string
	Attempt int // 第几次重试，从 1 开始
	Status  int // 触发重试的响应码；网络错误时为 0
	Err     error
	Wait    time.Duration
}

// RetryPolicy 是底层重试的完整配置。Sleep 与 OnRetry 可注入，测试用假时钟即可
// 让整套退避链路变成确定的同步流程，无需真实等待。
type RetryPolicy struct {
	MaxRetries    int
	BackoffFactor time.Duration
	MaxBackoff    time.Duration
	StatusCodes   map[int]struct{}
	Sleep         func(ctx context.Context, d time.Duration) error
	OnRetry       func(RetryEvent)
}

// statusSet 把状态码列表转换成集合。
func statusSet(codes []int) map[int]struct{} {
	set := make(map[int]struct{}, len(codes))
	for _, code := range codes {
		set[code] = struct{}{}
	}
	return set
}

// sleepContext 是默认可被 ctx 打断的等待实现。
func sleepContext(ctx context.Context, d time.Duration) error {
	if d <= 0 {
		return ctx.Err()
	}
	timer := time.NewTimer(d)
	defer timer.Stop()
	select {
	case <-ctx.Done():
		return ctx.Err()
	case <-timer.C:
		return nil
	}
}

// retryLimitKey 是请求级重试次数覆盖的 context key。
//
// 上传类接口用 NoRetries() 关闭底层重试，把重试交给上传管线的退避，
// 避免两层退避叠加放大请求量。
type retryLimitKey struct{}

// WithRetries 返回携带请求级重试次数覆盖的 context，用于关闭或收窄底层重试。
// 非负数才有意义，负数等价于 0。
func WithRetries(ctx context.Context, retries int) context.Context {
	if retries < 0 {
		retries = 0
	}
	return context.WithValue(ctx, retryLimitKey{}, retries)
}

// retriesFrom 读取请求级重试次数覆盖，未设置时返回 false。
func retriesFrom(ctx context.Context) (int, bool) {
	value, ok := ctx.Value(retryLimitKey{}).(int)
	return value, ok
}

// retryTransport 是承载底层重试的 RoundTripper。
//
// 重试判据与 urllib3 对齐：
//   - 状态码命中 StatusCodes；
//   - 传输层错误（连接失败、读超时），但 context 取消不重试；
//   - 重试前按 Retry-After（若有）或指数退避等待，等待期间响应体必须关闭以复用连接。
type retryTransport struct {
	base   http.RoundTripper
	policy RetryPolicy
}

var _ http.RoundTripper = (*retryTransport)(nil)

// newRetryTransport 在 base 之外套一层重试。base 为 nil 时使用 http.DefaultTransport。
func newRetryTransport(base http.RoundTripper, policy RetryPolicy) *retryTransport {
	if base == nil {
		base = http.DefaultTransport
	}
	if policy.Sleep == nil {
		policy.Sleep = sleepContext
	}
	if policy.BackoffFactor <= 0 {
		policy.BackoffFactor = DefaultBackoffFactor
	}
	if policy.MaxBackoff <= 0 {
		policy.MaxBackoff = DefaultMaxBackoff
	}
	if policy.StatusCodes == nil {
		policy.StatusCodes = statusSet(DefaultRetryStatusCodes)
	}
	if policy.MaxRetries < 0 {
		policy.MaxRetries = 0
	}
	return &retryTransport{base: base, policy: policy}
}

// RoundTrip 实现可重放的请求发送。请求体必须是可重放的（由 http.NewRequest 为
// bytes.Reader 自动填好 GetBody），否则第二次投递会直接失败，而不是发送半个体。
func (t *retryTransport) RoundTrip(req *http.Request) (*http.Response, error) {
	limit := t.policy.MaxRetries
	if override, ok := retriesFrom(req.Context()); ok {
		limit = override
	}

	for attempt := 0; ; attempt++ {
		attemptReq, err := replay(req, attempt)
		if err != nil {
			return nil, err
		}
		resp, err := t.base.RoundTrip(attemptReq)
		if err != nil {
			if !retryableTransportError(req.Context(), err) || attempt >= limit {
				return nil, err
			}
			wait := t.backoff(attempt+1, "")
			t.notify(RetryEvent{Method: req.Method, URL: req.URL.Redacted(), Attempt: attempt + 1, Err: err, Wait: wait})
			if err := t.policy.Sleep(req.Context(), wait); err != nil {
				return nil, err
			}
			continue
		}
		if _, retryable := t.policy.StatusCodes[resp.StatusCode]; !retryable || attempt >= limit {
			return resp, nil
		}
		wait := t.backoff(attempt+1, resp.Header.Get("Retry-After"))
		t.notify(RetryEvent{
			Method:  req.Method,
			URL:     req.URL.Redacted(),
			Attempt: attempt + 1,
			Status:  resp.StatusCode,
			Wait:    wait,
		})
		// 必须排空并关闭响应体，否则连接无法回到池里复用。
		drainAndClose(resp.Body)
		if err := t.policy.Sleep(req.Context(), wait); err != nil {
			return nil, err
		}
	}
}

// CloseIdleConnections 转发给底层连接池，供 Factory.Close 释放空闲连接。
func (t *retryTransport) CloseIdleConnections() {
	if closer, ok := t.base.(interface{ CloseIdleConnections() }); ok {
		closer.CloseIdleConnections()
	}
}

// replay 为第 attempt 次发送准备请求副本。首次直接克隆原请求，之后通过 GetBody
// 重放请求体；缺少 GetBody 时明确失败，避免发出残缺请求。
func replay(req *http.Request, attempt int) (*http.Request, error) {
	clone := req.Clone(req.Context())
	if attempt == 0 || req.Body == nil || req.Body == http.NoBody {
		return clone, nil
	}
	if req.GetBody == nil {
		return nil, fmt.Errorf("%w: request body is not replayable", ErrInvalidResponse)
	}
	body, err := req.GetBody()
	if err != nil {
		return nil, err
	}
	clone.Body = body
	return clone, nil
}

// backoff 计算第 attempt 次重试前的等待时长（attempt 从 1 开始）。
//
// 优先使用服务端给出的 Retry-After；否则按 factor * 2^(attempt-1) 指数增长，
// 两者都受 MaxBackoff 约束，避免异常响应头导致长时间挂起。
func (t *retryTransport) backoff(attempt int, retryAfter string) time.Duration {
	if hint := parseRetryAfter(retryAfter, time.Now()); hint > 0 {
		return min(hint, t.policy.MaxBackoff)
	}
	wait := t.policy.BackoffFactor
	for i := 1; i < attempt; i++ {
		wait *= 2
		if wait >= t.policy.MaxBackoff {
			return t.policy.MaxBackoff
		}
	}
	return wait
}

// notify 把重试事件交给可选观测回调。
func (t *retryTransport) notify(event RetryEvent) {
	if t.policy.OnRetry != nil {
		t.policy.OnRetry(event)
	}
}

// retryableTransportError 判断传输层错误是否值得重放。
//
// 调用方已终止（ctx 取消或 deadline）时不重试；其余网络错误（连接失败、读写超时）
// 均保留重试，与 Python urllib3 对 connect/read 错误的重试语义一致。
func retryableTransportError(ctx context.Context, err error) bool {
	if err == nil {
		return false
	}
	if errors.Is(err, context.Canceled) {
		return false
	}
	if ctx.Err() != nil {
		return false
	}
	return true
}

// drainAndClose 丢弃响应体并关闭，用于重试前回收连接。
func drainAndClose(body io.ReadCloser) {
	if body == nil {
		return
	}
	_, _ = io.Copy(io.Discard, io.LimitReader(body, 64<<10))
	_ = body.Close()
}
