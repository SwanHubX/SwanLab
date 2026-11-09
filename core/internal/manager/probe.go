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
	// controlProbePath 是控制面存活端点；Server#1168 上线后可切到 /api/health。
	controlProbePath = "/api/"
	// dataProbePath 是数据面存活端点，返回静态 200。
	dataProbePath = "/api/house/health"
	// probeTimeout 是单次探测预算；启动预检与周期探测共用。
	probeTimeout = 2 * time.Second
)

// probeClient 与认证客户端隔离：不携带 API Key，也不参与认证重试。
var probeClient = &http.Client{Timeout: probeTimeout}

// siteProbe 把站点地址与路径包装成可注入的探针函数。
func siteProbe(apiHost, path string) health.ProbeFunc {
	url := strings.TrimRight(apiHost, "/") + path
	return func(ctx context.Context) bool { return probeReachable(ctx, url) }
}

// probeNow 用独立的 2s 预算执行一次探测，供启动预检使用。
func probeNow(ctx context.Context, apiHost, path string) bool {
	ctx, cancel := context.WithTimeout(ctx, probeTimeout)
	defer cancel()
	return siteProbe(apiHost, path)(ctx)
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
	// 消费 body 以便复用连接；探测响应很小，限量即可。
	_, _ = io.Copy(io.Discard, io.LimitReader(resp.Body, 4<<10))
	_ = resp.Body.Close()
	return resp.StatusCode != http.StatusTooManyRequests && resp.StatusCode < 500
}

// probePath 返回配置的探针路径；空串使用默认路径。
func probePath(configured, fallback string) string {
	if configured == "" {
		return fallback
	}
	return configured
}

// secondsDuration 把协议中的秒数转换为 Duration；<=0 保持 0，由 health 使用默认值。
func secondsDuration(seconds float64) time.Duration {
	if seconds <= 0 {
		return 0
	}
	return time.Duration(seconds * float64(time.Second))
}

// probeStartup 依次探测控制面与数据面；控制面不可达时不再探测数据面。
func probeStartup(ctx context.Context, apiHost, controlPath, dataPath string) (control, data bool) {
	if !probeNow(ctx, apiHost, controlPath) {
		return false, false
	}
	return true, probeNow(ctx, apiHost, dataPath)
}

// siteKey 用于在健康注册表中按 api_host + api_key 摘要共享探针，
// 避免把密钥本身放进内存键或日志。
func siteKey(apiHost, apiKey string) string {
	sum := sha256.Sum256([]byte(apiHost + "\x00" + apiKey))
	return hex.EncodeToString(sum[:16])
}
