package client

import (
	"encoding/json"
	"errors"
	"fmt"
	"net/http"
	"strconv"
	"strings"
	"time"
)

// traceHeader 是服务端回传的链路追踪头，写入 APIError 供排查。
const traceHeader = "traceid"

// 领域哨兵错误，便于调用方用 errors.Is 分流而不必解析字符串。
var (
	// ErrNotAuthenticated 表示凭证被服务端拒绝（401/403）。
	ErrNotAuthenticated = errors.New("client: authentication failed")
	// ErrInvalidResponse 表示响应体不符合协议（字段缺失、类型错误、数量不匹配）。
	// 这类错误重试不会恢复，调用方应直接终结该批记录。
	ErrInvalidResponse = errors.New("client: invalid response payload")
	// ErrHouseRejected 表示 House 接口未以 200 OK 确认写入（典型形态是
	// 202 Accepted + "session is invalid"）。后端没有落库，重试同一会话也不会成功。
	ErrHouseRejected = errors.New("client: house rejected the write")
	// ErrMissingAPIKey 表示未提供 api_key。
	ErrMissingAPIKey = errors.New("client: api_key must not be empty")
	// ErrMissingAPIHost 表示未提供 api_host。
	ErrMissingAPIHost = errors.New("client: api_host must not be empty")
)

// APIError 是与云端交互失败的统一错误类型。
//
// 字段与 Python exceptions.ApiError 对齐：业务错误码、人类可读消息、
// 链路 TraceID，以及显式声明的重试语义。Retryable 的判定规则：
//   - 显式指定优先（协议错误恒为 false，与 HTTP 状态码无关）；
//   - 否则 429 与 5xx 可重试，其余 4xx 视为确定性失败；
//   - 没有响应（网络错误、超时）保留重试。
type APIError struct {
	Method     string
	URL        string
	StatusCode int // 0 表示请求未拿到响应（连接失败、超时、协议错误）
	Code       string
	Message    string
	TraceID    string
	Retryable  bool
	// RetryAfter 是服务端显式要求的冷却窗口（Retry-After 头），0 表示未提供。
	RetryAfter time.Duration
	// Err 保留底层原因（网络错误、协议错误哨兵），支持 errors.Is 判定。
	Err error
}

// Error 对齐 Python 的 "[code] message | Trace ID | METHOD url" 排版，便于对照排查。
func (e *APIError) Error() string {
	var b strings.Builder
	b.WriteString("API Request Failed: [")
	b.WriteString(e.Code)
	b.WriteString("] ")
	b.WriteString(e.Message)
	if e.TraceID != "" {
		b.WriteString(" | Trace ID: ")
		b.WriteString(e.TraceID)
	}
	if e.Method != "" || e.URL != "" {
		b.WriteString(" | ")
		b.WriteString(e.Method)
		b.WriteString(" ")
		b.WriteString(e.URL)
	}
	return b.String()
}

// Unwrap 暴露底层原因，支持 errors.Is(err, ErrNotAuthenticated) 等判定。
func (e *APIError) Unwrap() error { return e.Err }

// IsRetryable 判断错误是否属于「保留记录稍后重试」的一类。
// 非 APIError（例如 context 取消）一律视为不可重试，避免无意义的重放。
func IsRetryable(err error) bool {
	var apiErr *APIError
	if errors.As(err, &apiErr) {
		return apiErr.Retryable
	}
	return false
}

// IsAuthenticationError 判断错误是否由凭证被拒绝导致。
func IsAuthenticationError(err error) bool {
	if errors.Is(err, ErrNotAuthenticated) {
		return true
	}
	var apiErr *APIError
	return errors.As(err, &apiErr) && (apiErr.StatusCode == http.StatusUnauthorized || apiErr.StatusCode == http.StatusForbidden)
}

// IsInvalidResponse 判断错误是否属于响应不符合协议。
func IsInvalidResponse(err error) bool { return errors.Is(err, ErrInvalidResponse) }

// IsHouseRejected 判断错误是否属于「House 未以 200 确认写入」。
// 这是会话级失败，上层应显式告警并终结该批记录，而不是继续退避重试。
func IsHouseRejected(err error) bool { return errors.Is(err, ErrHouseRejected) }

// IsNotFound 判断错误是否为 404。
func IsNotFound(err error) bool { return StatusCode(err) == http.StatusNotFound }

// StatusCode 返回 APIError 的 HTTP 状态码，非 APIError 返回 0。
func StatusCode(err error) int {
	var apiErr *APIError
	if errors.As(err, &apiErr) {
		return apiErr.StatusCode
	}
	return 0
}

// newTransportError 把底层网络错误包装成「可重试」的 APIError。
// 调用方此时已完成底层重试，仍需保留重试语义：管线层会退避后重新投递。
func newTransportError(method, url string, err error) *APIError {
	return &APIError{
		Method:    method,
		URL:       url,
		Code:      "network_error",
		Message:   err.Error(),
		Retryable: true,
		Err:       err,
	}
}

// newAPIError 依据 HTTP 响应构造 APIError，业务码与消息从响应体解析，解析不出时
// 回落到状态码与状态文本（对齐 Python decode_error_response）。
func newAPIError(method, endpoint string, status int, header http.Header, body []byte) *APIError {
	statusText := http.StatusText(status)
	code, message := decodeErrorBody(body, status, statusText)
	traceID := header.Get(traceHeader)
	apiErr := &APIError{
		Method:     method,
		URL:        endpoint,
		StatusCode: status,
		Code:       code,
		Message:    message,
		TraceID:    traceID,
		Retryable:  retryableStatus(status),
		RetryAfter: parseRetryAfter(header.Get("Retry-After"), time.Now()),
	}
	if status == http.StatusUnauthorized || status == http.StatusForbidden {
		apiErr.Err = ErrNotAuthenticated
	}
	return apiErr
}

// NewProtocolError 构造不可重试的协议错误：响应无法按契约解读。
// 供 internal/api 在字段缺失、数量不匹配等场景复用，保证错误分类一致。
func NewProtocolError(method, url, code, message string, cause error) *APIError {
	return &APIError{
		Method:    method,
		URL:       url,
		Code:      code,
		Message:   message,
		Retryable: false,
		Err:       errors.Join(ErrInvalidResponse, cause),
	}
}

// HouseRejected 把「House 未以 200 确认写入」的响应转成不可重试的业务错误。
//
// 典型形态是 202 Accepted + "session is invalid"/"FlagId is invalid"：后端没有写库，
// 重试同一会话也不会成功，必须显式判定失败并让上层终结该批记录。
func (r *Response) HouseRejected(code string) *APIError {
	message := strings.TrimSpace(string(r.Body))
	if message == "" {
		message = fmt.Sprintf("house returned status %d instead of 200", r.Status)
	}
	return &APIError{
		Method:     r.Method,
		URL:        r.Endpoint,
		StatusCode: r.Status,
		Code:       code,
		Message:    truncate(message, 512),
		TraceID:    r.Header.Get(traceHeader),
		Retryable:  false,
		Err:        ErrHouseRejected,
	}
}

// decodeErrorBody 解析后端错误响应体中的 {code, message}；不是 JSON 对象时回落到
// 状态码与状态文本，保证错误消息始终可读。
func decodeErrorBody(body []byte, status int, statusText string) (string, string) {
	code := strconv.Itoa(status)
	message := statusText
	if len(body) == 0 {
		return code, message
	}
	var payload map[string]any
	if err := json.Unmarshal(body, &payload); err != nil || payload == nil {
		// 非 JSON 错误体：把原始文本作为消息，便于排查网关注入的 HTML 错误页。
		if text := strings.TrimSpace(string(body)); text != "" {
			message = truncate(text, 512)
		}
		return code, message
	}
	if v, ok := payload["code"]; ok && v != nil {
		code = fmt.Sprint(v)
	}
	if v, ok := payload["message"]; ok && v != nil {
		message = fmt.Sprint(v)
	}
	return code, message
}

// retryableStatus 实现「429 与 5xx 可重试、其余 4xx 确定性失败」的分类。
func retryableStatus(status int) bool {
	return status == http.StatusTooManyRequests || status >= http.StatusInternalServerError
}

// parseRetryAfter 解析 Retry-After（秒数或 HTTP-date），无法解析时返回 0。
func parseRetryAfter(value string, now time.Time) time.Duration {
	value = strings.TrimSpace(value)
	if value == "" {
		return 0
	}
	if seconds, err := strconv.Atoi(value); err == nil {
		if seconds <= 0 {
			return 0
		}
		return time.Duration(seconds) * time.Second
	}
	if at, err := http.ParseTime(value); err == nil {
		if wait := at.Sub(now); wait > 0 {
			return wait
		}
	}
	return 0
}

// truncate 截断过长文本，避免把整页 HTML 塞进日志。
func truncate(text string, limit int) string {
	if len(text) <= limit {
		return text
	}
	return text[:limit] + "..."
}
