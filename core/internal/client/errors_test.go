package client

import (
	"context"
	"errors"
	"io"
	"net/http"
	"strings"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
)

// TestHouseRejectedError 覆盖「House 未以 200 确认写入」的错误构造：
// 保留后端原文与 TraceID，且显式声明不可重试。
func TestHouseRejectedError(t *testing.T) {
	header := http.Header{}
	header.Set(traceHeader, "trace-9")
	res := &Response{
		Method:   http.MethodPost,
		Endpoint: "https://api.example.com/api/house/metrics",
		Status:   http.StatusAccepted,
		Header:   header,
		Body:     []byte("session is invalid"),
	}

	err := res.HouseRejected("session_rejected")
	if err.StatusCode != http.StatusAccepted || err.Code != "session_rejected" {
		t.Fatalf("error = %+v, want 202 session_rejected", err)
	}
	if err.Retryable || !IsHouseRejected(err) {
		t.Fatalf("error = %+v, want non-retryable house rejection", err)
	}
	if err.TraceID != "trace-9" || err.Message != "session is invalid" {
		t.Fatalf("error = %+v, want trace id and backend body preserved", err)
	}

	empty := (&Response{Status: http.StatusAccepted, Header: http.Header{}}).HouseRejected("heartbeat_rejected")
	if !strings.Contains(empty.Message, "instead of 200") {
		t.Fatalf("empty body message = %q, want explicit status hint", empty.Message)
	}
}

// TestParseRetryAfter 覆盖 Retry-After 的两种合法形态与非法输入。
func TestParseRetryAfter(t *testing.T) {
	now := clienttest.FixedNow
	future := now.Add(30 * time.Second).Format(http.TimeFormat)
	past := now.Add(-time.Minute).Format(http.TimeFormat)

	cases := []struct {
		in   string
		want time.Duration
	}{
		{in: "", want: 0},
		{in: "  ", want: 0},
		{in: "5", want: 5 * time.Second},
		{in: "0", want: 0},
		{in: "-3", want: 0},
		{in: past, want: 0},
		{in: "not-a-date", want: 0},
	}
	for _, tc := range cases {
		if got := parseRetryAfter(tc.in, now); got != tc.want {
			t.Fatalf("parseRetryAfter(%q) = %v, want %v", tc.in, got, tc.want)
		}
	}

	// HTTP-date 只有秒级精度，允许 1s 误差。
	got := parseRetryAfter(future, now)
	if got < 29*time.Second || got > 30*time.Second {
		t.Fatalf("parseRetryAfter(http-date) = %v, want ~30s", got)
	}
}

// TestDecodeErrorBody 覆盖错误响应体解析：JSON 取值、非 JSON 回落、截断。
func TestDecodeErrorBody(t *testing.T) {
	code, message := decodeErrorBody([]byte(`{"code": "expired", "message": "sid expired"}`), 401, "Unauthorized")
	if code != "expired" || message != "sid expired" {
		t.Fatalf("decoded = (%q, %q), want (expired, sid expired)", code, message)
	}

	code, message = decodeErrorBody([]byte("<html>gateway timeout</html>"), http.StatusGatewayTimeout, "Gateway Timeout")
	if code != "504" || !strings.Contains(message, "gateway timeout") {
		t.Fatalf("fallback = (%q, %q), want status code and raw text", code, message)
	}

	code, message = decodeErrorBody(nil, http.StatusBadRequest, "Bad Request")
	if code != "400" || message != "Bad Request" {
		t.Fatalf("empty body = (%q, %q), want status code and text", code, message)
	}

	_, truncated := decodeErrorBody([]byte(strings.Repeat("x", 1024)), 500, "Internal Server Error")
	if len(truncated) > 520 || !strings.HasSuffix(truncated, "...") {
		t.Fatalf("long body must be truncated, got len=%d", len(truncated))
	}
}

// TestAPIErrorFormatting 覆盖错误文案排版与 unwrap 链。
func TestAPIErrorFormatting(t *testing.T) {
	apiErr := &APIError{
		Method:     http.MethodGet,
		URL:        "https://api.example.com/api/x",
		StatusCode: http.StatusUnauthorized,
		Code:       "unauthorized",
		Message:    "invalid api key",
		TraceID:    "trace-1",
		Err:        ErrNotAuthenticated,
	}
	text := apiErr.Error()
	for _, want := range []string{"[unauthorized]", "invalid api key", "Trace ID: trace-1", "GET https://api.example.com/api/x"} {
		if !strings.Contains(text, want) {
			t.Fatalf("error text %q missing %q", text, want)
		}
	}
	if !IsAuthenticationError(apiErr) {
		t.Fatal("401 APIError must be classified as authentication error")
	}
}

// TestRetryableClassification 覆盖重试语义分类表：429/5xx 可重试，其余 4xx 确定性失败。
func TestRetryableClassification(t *testing.T) {
	cases := []struct {
		name   string
		status int
		want   bool
	}{
		{name: "bad request", status: http.StatusBadRequest, want: false},
		{name: "unauthorized", status: http.StatusUnauthorized, want: false},
		{name: "not found", status: http.StatusNotFound, want: false},
		{name: "rate limited", status: http.StatusTooManyRequests, want: true},
		{name: "internal", status: http.StatusInternalServerError, want: true},
		{name: "unavailable", status: http.StatusServiceUnavailable, want: true},
	}
	for _, tc := range cases {
		err := &APIError{StatusCode: tc.status, Retryable: retryableStatus(tc.status)}
		if got := IsRetryable(err); got != tc.want {
			t.Fatalf("%s: IsRetryable = %v, want %v", tc.name, got, tc.want)
		}
	}
	// 协议错误与鉴权错误显式声明不可重试，与状态码无关。
	if IsRetryable(NewProtocolError("GET", "/x", "bad", "bad", nil)) {
		t.Fatal("protocol error must not be retryable")
	}
	if IsRetryable(errors.New("plain")) {
		t.Fatal("non-APIError must not be retryable")
	}
	if IsRetryable(nil) {
		t.Fatal("nil must not be retryable")
	}
}

// TestRetryableTransportError 覆盖传输层错误的重试判据。
func TestRetryableTransportError(t *testing.T) {
	background := context.Background()
	if retryableTransportError(background, nil) {
		t.Fatal("nil error must not be retryable")
	}
	if retryableTransportError(background, context.Canceled) {
		t.Fatal("canceled context must not be retried")
	}
	if !retryableTransportError(background, errors.New("connection reset by peer")) {
		t.Fatal("network error must be retryable")
	}
	// 客户端超时（http.Client.Timeout）不是调用方主动终止，应与 Python read timeout 一样重试。
	if !retryableTransportError(background, context.DeadlineExceeded) {
		t.Fatal("client-side deadline must be retryable")
	}
	// 调用方 ctx 已终止时，任何错误都不再重试。
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if retryableTransportError(ctx, errors.New("connection reset by peer")) {
		t.Fatal("transport error after caller cancellation must not be retried")
	}
}

// TestSleepContext 覆盖默认等待实现可被 context 打断。
func TestSleepContext(t *testing.T) {
	if err := sleepContext(context.Background(), 5*time.Millisecond); err != nil {
		t.Fatalf("sleep: %v", err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if err := sleepContext(ctx, time.Minute); !errors.Is(err, context.Canceled) {
		t.Fatalf("sleep = %v, want context.Canceled", err)
	}
}

// TestReplayRejectsNonReplayableBody 覆盖无法重放的请求体：明确失败而不是发出残缺请求。
func TestReplayRejectsNonReplayableBody(t *testing.T) {
	req, err := http.NewRequest(http.MethodPost, "https://api.example.com/api/x", io.NopCloser(strings.NewReader("body")))
	if err != nil {
		t.Fatalf("new request: %v", err)
	}
	req.GetBody = nil
	if _, err := replay(req, 0); err != nil {
		t.Fatalf("first attempt must not need GetBody: %v", err)
	}
	if _, err := replay(req, 1); !errors.Is(err, ErrInvalidResponse) {
		t.Fatalf("replay = %v, want ErrInvalidResponse", err)
	}
}

// TestProfileValidate 覆盖鉴权响应字段校验。
func TestProfileValidate(t *testing.T) {
	valid := Profile{UID: 1, Username: "alice", CreatedAt: "2026-01-01T00:00:00Z"}
	if err := valid.validate(http.MethodGet, "/api/auth/verify"); err != nil {
		t.Fatalf("valid profile rejected: %v", err)
	}
	cases := []struct {
		name    string
		profile Profile
	}{
		{name: "zero uid", profile: Profile{Username: "alice", CreatedAt: "x"}},
		{name: "empty username", profile: Profile{UID: 1, CreatedAt: "x"}},
		{name: "empty created at", profile: Profile{UID: 1, Username: "alice"}},
	}
	for _, tc := range cases {
		if err := tc.profile.validate("GET", "/api/auth/verify"); !IsInvalidResponse(err) {
			t.Fatalf("%s: error = %v, want invalid response", tc.name, err)
		}
	}
}

// TestParseExpiredAt 覆盖登录响应过期时间的多种排版。
func TestParseExpiredAt(t *testing.T) {
	want := time.Date(2999, 1, 1, 0, 0, 0, 0, time.UTC)
	cases := []string{
		"2999-01-01T00:00:00.000000Z",
		"2999-01-01T00:00:00Z",
		"2999-01-01T00:00:00",
	}
	for _, tc := range cases {
		got, err := parseExpiredAt(tc)
		if err != nil {
			t.Fatalf("parseExpiredAt(%q): %v", tc, err)
		}
		if !got.Equal(want) {
			t.Fatalf("parseExpiredAt(%q) = %v, want %v", tc, got, want)
		}
	}
	if _, err := parseExpiredAt(""); err == nil {
		t.Fatal("empty expiredAt must be rejected")
	}
	if _, err := parseExpiredAt("tomorrow"); err == nil {
		t.Fatal("unparsable expiredAt must be rejected")
	}
}

// TestEncodedBodyPassthrough 覆盖 []byte 请求体直通（供预编码 payload 复用）。
func TestEncodedBodyPassthrough(t *testing.T) {
	raw := []byte(`{"a":1}`)
	encoded, err := encodeBody(raw)
	if err != nil {
		t.Fatalf("encode body: %v", err)
	}
	if string(encoded) != string(raw) {
		t.Fatalf("encoded = %s, want passthrough", encoded)
	}
	if encoded, err := encodeBody(nil); err != nil || encoded != nil {
		t.Fatalf("nil body = (%v, %v), want (nil, nil)", encoded, err)
	}
	if _, err := encodeBody(func() {}); err == nil {
		t.Fatal("unencodable body must be rejected")
	}
}
