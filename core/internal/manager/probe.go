// 本文件提供站点探针的 HTTP 实现与站点标识；探针不经过认证，也不使用认证客户端。
package manager

import (
	"context"
	"crypto/sha256"
	"encoding/hex"
	"io"
	"net/http"
	"strings"
	"time"

	"github.com/swanhubx/swanlab/core/internal/cron/health"
)

const (
	// serverProbePath 是 swanlab-server 存活端点，返回静态 200。
	serverProbePath = "/api/health"
	// houseProbePath 是 swanlab-house 存活端点，返回静态 200。
	houseProbePath = "/api/house/health"
	// probeTimeout 是启动预检的单次探测预算；周期探测使用 health_probe_timeout。
	probeTimeout = 2 * time.Second
)

// probeClient 与认证客户端隔离：不携带 API Key，也不参与认证重试。
// 不设 Client.Timeout：预算统一由请求 context 控制，避免截断 health_probe_timeout。
var probeClient = &http.Client{}

// probeURL 拼接 api_host 与探针路径。
func probeURL(apiHost, path string) string {
	return strings.TrimRight(apiHost, "/") + path
}

// siteProbe 把站点地址与路径包装成可注入的探针函数。
func siteProbe(apiHost, path string) health.ProbeFunc {
	url := probeURL(apiHost, path)
	return func(ctx context.Context) bool { return probeReachable(ctx, url) }
}

// probeStartup 用固定预算依次探测 swanlab-server 与 swanlab-house；server 不可达时不再探测 house。
func probeStartup(ctx context.Context, server, house health.ProbeFunc) (serverOK, houseOK bool) {
	probe := func(fn health.ProbeFunc) bool {
		pctx, cancel := context.WithTimeout(ctx, probeTimeout)
		defer cancel()
		return fn(pctx)
	}
	if !probe(server) {
		return false, false
	}
	return true, probe(house)
}

// probeReachable 只判断服务是否可达：2xx/3xx/4xx 视为可达（4xx 是确定性拒绝，
// 不能驱动暂停循环），transport error、5xx 与 429 视为不可达。
func probeReachable(ctx context.Context, url string) bool {
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, url, nil)
	if err != nil {
		return false
	}
	resp, err := probeClient.Do(req)
	if err != nil {
		return false
	}
	// 读尽 body 以便复用连接；读取受请求 ctx 的预算约束。
	_, _ = io.Copy(io.Discard, resp.Body)
	_ = resp.Body.Close()
	return resp.StatusCode != http.StatusTooManyRequests && resp.StatusCode < 500
}

// secondsDuration 把协议中的秒数转换为 Duration；<=0 保持 0，由 health 使用默认值。
func secondsDuration(seconds float64) time.Duration {
	if seconds <= 0 {
		return 0
	}
	return time.Duration(seconds * float64(time.Second))
}

// siteKey 按 api_host + api_key 摘要生成站点标识，密钥本身不进入内存键或日志。
func siteKey(apiHost, apiKey string) string {
	sum := sha256.Sum256([]byte(apiHost + "\x00" + apiKey))
	return hex.EncodeToString(sum[:16])
}
