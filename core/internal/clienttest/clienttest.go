// Package clienttest 提供 internal/client 与 internal/api 共用的测试辅助：
// httptest 假后端与可注入的假时钟。只被 *_test.go 引用，不进入生产路径。
package clienttest

import (
	"context"
	"encoding/json"
	"io"
	"net/http"
	"net/http/httptest"
	"net/url"
	"sync"
	"testing"
	"time"
)

// FixedNow 是测试用的固定时间，带微秒以便断言 createdAt / finishedAt 的 6 位小数排版。
var FixedNow = time.Date(2026, 3, 10, 19, 2, 3, 123400000, time.UTC)

// Request 是被 Cloud 记录下来的请求快照。
type Request struct {
	Method string
	Path   string
	Query  url.Values
	Header http.Header
	Body   []byte
}

// Cloud 是运行在 httptest 上的极简 SwanLab 后端，按 "METHOD path" 注册路由。
//
// 所有请求都会被记录，测试据此断言请求头（鉴权、SDK 版本）、请求体与重试次数；
// 未注册的路径返回 404，用于覆盖「旧后端缺端点」等分支。
type Cloud struct {
	t        *testing.T
	server   *httptest.Server
	mu       sync.Mutex
	routes   map[string]http.HandlerFunc
	requests []Request
	close    sync.Once
}

// NewCloud 启动一个本地假后端，测试结束时自动关闭。
func NewCloud(t *testing.T) *Cloud {
	t.Helper()
	cloud := &Cloud{t: t, routes: make(map[string]http.HandlerFunc)}
	cloud.server = httptest.NewServer(http.HandlerFunc(cloud.serve))
	t.Cleanup(cloud.Close)
	return cloud
}

// URL 返回假后端的根地址。
func (c *Cloud) URL() string { return c.server.URL }

// Handle 注册路由处理函数。
func (c *Cloud) Handle(method, path string, handler http.HandlerFunc) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.routes[method+" "+path] = handler
}

// ServeJSON 注册返回固定状态码与 JSON 体的路由。payload 为 nil 时只写状态码。
func (c *Cloud) ServeJSON(method, path string, status int, payload any) {
	c.Handle(method, path, func(w http.ResponseWriter, _ *http.Request) {
		WriteJSON(w, status, payload)
	})
}

// Close 关闭假后端。可重复调用，便于模拟「后端不可达」后仍能被 Cleanup 收尾。
func (c *Cloud) Close() {
	c.close.Do(c.server.Close)
}

// Requests 返回按顺序记录的请求快照。
func (c *Cloud) Requests() []Request {
	c.mu.Lock()
	defer c.mu.Unlock()
	return append([]Request(nil), c.requests...)
}

// RequestsFor 返回命中指定 METHOD+path 的请求快照。
func (c *Cloud) RequestsFor(method, path string) []Request {
	var matched []Request
	for _, req := range c.Requests() {
		if req.Method == method && req.Path == path {
			matched = append(matched, req)
		}
	}
	return matched
}

// CountFor 返回命中指定 METHOD+path 的请求次数。
func (c *Cloud) CountFor(method, path string) int { return len(c.RequestsFor(method, path)) }

// serve 记录请求后分发到路由表。
func (c *Cloud) serve(w http.ResponseWriter, r *http.Request) {
	body, _ := io.ReadAll(r.Body)
	c.mu.Lock()
	c.requests = append(c.requests, Request{
		Method: r.Method,
		Path:   r.URL.Path,
		Query:  r.URL.Query(),
		Header: r.Header.Clone(),
		Body:   body,
	})
	handler := c.routes[r.Method+" "+r.URL.Path]
	c.mu.Unlock()

	if handler == nil {
		WriteJSON(w, http.StatusNotFound, map[string]string{"code": "not_found", "message": "no route " + r.URL.Path})
		return
	}
	handler(w, r)
}

// WriteJSON 输出 JSON 响应。
func WriteJSON(w http.ResponseWriter, status int, payload any) {
	w.Header().Set("Content-Type", "application/json")
	w.WriteHeader(status)
	if payload == nil {
		return
	}
	_ = json.NewEncoder(w).Encode(payload)
}

// Clock 是可注入的时钟：Now 固定在 FixedNow，Sleep 只记录时长不真实等待，
// 使退避时间线可以被精确断言。
type Clock struct {
	mu     sync.Mutex
	now    time.Time
	sleeps []time.Duration
}

// NewClock 返回从 FixedNow 起步的假时钟。
func NewClock() *Clock { return &Clock{now: FixedNow} }

// Now 返回当前假时间。
func (c *Clock) Now() time.Time {
	c.mu.Lock()
	defer c.mu.Unlock()
	return c.now
}

// Sleep 记录等待时长并推进假时间。
func (c *Clock) Sleep(_ context.Context, d time.Duration) error {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.sleeps = append(c.sleeps, d)
	c.now = c.now.Add(d)
	return nil
}

// Sleeps 返回累计的等待时长序列。
func (c *Clock) Sleeps() []time.Duration {
	c.mu.Lock()
	defer c.mu.Unlock()
	return append([]time.Duration(nil), c.sleeps...)
}
