package client

import (
	"errors"
	"io"
	"net/http"
	"net/http/httptest"
	"sync"
	"sync/atomic"
	"testing"
	"time"
)

func newHTTP2TestServer(t *testing.T, handler http.Handler) *httptest.Server {
	t.Helper()
	server := httptest.NewUnstartedServer(handler)
	server.EnableHTTP2 = true
	server.StartTLS()
	t.Cleanup(server.Close)
	return server
}

func cloneHTTP2TestTransport(t *testing.T, server *httptest.Server) *http.Transport {
	t.Helper()
	base, ok := server.Client().Transport.(*http.Transport)
	if !ok {
		t.Fatalf("unexpected server client transport type: %T", server.Client().Transport)
	}
	return base.Clone()
}

// TestTransportPoolReusesFirstShardWhenSequential 验证串行流量始终复用分片 0 的已建 HTTP/2 连接，
// 不为低并发请求触发额外 TLS 握手。
func TestTransportPoolReusesFirstShardWhenSequential(t *testing.T) {
	const (
		shardCount = 4
		requests   = 5
		http2Major = 2
	)
	var mu sync.Mutex
	remotes := make([]string, 0, requests)
	protos := make([]int, 0, requests)

	server := newHTTP2TestServer(t, http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		mu.Lock()
		remotes = append(remotes, r.RemoteAddr)
		protos = append(protos, r.ProtoMajor)
		mu.Unlock()
		w.WriteHeader(http.StatusOK)
	}))

	pool := newTransportPool(shardCount, func() http.RoundTripper {
		return cloneHTTP2TestTransport(t, server)
	})
	t.Cleanup(pool.CloseIdleConnections)

	httpClient := &http.Client{Transport: pool}
	for i := 0; i < requests; i++ {
		resp, err := httpClient.Get(server.URL)
		if err != nil {
			t.Fatalf("request %d failed: %v", i, err)
		}
		_, _ = io.Copy(io.Discard, resp.Body)
		if closeErr := resp.Body.Close(); closeErr != nil {
			t.Fatalf("close body %d: %v", i, closeErr)
		}
	}

	if len(protos) != requests {
		t.Fatalf("recorded requests = %d, want %d", len(protos), requests)
	}
	for i, proto := range protos {
		if proto != http2Major {
			t.Fatalf("request %d ProtoMajor = %d, want %d", i, proto, http2Major)
		}
	}
	for i := 1; i < len(remotes); i++ {
		if remotes[i] != remotes[0] {
			t.Fatalf("request %d used remote %s, want reuse of %s", i, remotes[i], remotes[0])
		}
	}
}

// TestTransportPoolSpreadsConcurrentHTTP2AcrossConnections 验证并发 HTTP/2 请求
// 能够横向展开到分片数对应的多条独立 TCP 连接。
func TestTransportPoolSpreadsConcurrentHTTP2AcrossConnections(t *testing.T) {
	const (
		shardCount  = 4
		http2Major  = 2
		waitTimeout = 5 * time.Second
	)
	var arrived atomic.Int32
	allArrived := make(chan struct{})
	release := make(chan struct{})
	t.Cleanup(func() {
		select {
		case <-release:
		default:
			close(release)
		}
	})

	var mu sync.Mutex
	uniqueRemotes := make(map[string]struct{}, shardCount)
	protos := make([]int, 0, shardCount)

	server := newHTTP2TestServer(t, http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		mu.Lock()
		uniqueRemotes[r.RemoteAddr] = struct{}{}
		protos = append(protos, r.ProtoMajor)
		mu.Unlock()

		if int(arrived.Add(1)) == shardCount {
			close(allArrived)
		}
		select {
		case <-release:
		case <-r.Context().Done():
		}
		w.WriteHeader(http.StatusOK)
	}))

	pool := newTransportPool(shardCount, func() http.RoundTripper {
		return cloneHTTP2TestTransport(t, server)
	})
	t.Cleanup(pool.CloseIdleConnections)

	httpClient := &http.Client{Transport: pool}
	errs := make([]error, shardCount)
	var wg sync.WaitGroup
	for i := 0; i < shardCount; i++ {
		wg.Add(1)
		go func(idx int) {
			defer wg.Done()
			resp, err := httpClient.Get(server.URL)
			if err != nil {
				errs[idx] = err
				return
			}
			_, _ = io.Copy(io.Discard, resp.Body)
			if closeErr := resp.Body.Close(); closeErr != nil {
				errs[idx] = closeErr
			}
		}(i)
	}

	select {
	case <-allArrived:
		close(release)
	case <-time.After(waitTimeout):
		close(release)
		wg.Wait()
		t.Fatalf("timed out waiting for %d concurrent requests (arrived=%d)", shardCount, arrived.Load())
	}
	wg.Wait()

	for i, err := range errs {
		if err != nil {
			t.Fatalf("concurrent request %d failed: %v", i, err)
		}
	}
	if len(protos) != shardCount {
		t.Fatalf("recorded requests = %d, want %d", len(protos), shardCount)
	}
	for i, proto := range protos {
		if proto != http2Major {
			t.Fatalf("request %d ProtoMajor = %d, want %d", i, proto, http2Major)
		}
	}
	if got := len(uniqueRemotes); got != shardCount {
		t.Fatalf("unique TCP remote addrs = %d, want %d (remotes=%v)", got, shardCount, uniqueRemotes)
	}
}

type blockingBody struct {
	closed atomic.Bool
}

func (*blockingBody) Read([]byte) (int, error) { return 0, io.EOF }

func (b *blockingBody) Close() error {
	b.closed.Store(true)
	return nil
}

type barrierRoundTripper struct {
	maxConcurrent atomic.Int64
	active        atomic.Int64
	gate          <-chan struct{}
}

func (b *barrierRoundTripper) RoundTrip(req *http.Request) (*http.Response, error) {
	cur := b.active.Add(1)
	for {
		prev := b.maxConcurrent.Load()
		if cur <= prev || b.maxConcurrent.CompareAndSwap(prev, cur) {
			break
		}
	}
	<-b.gate
	b.active.Add(-1)
	return &http.Response{
		StatusCode: http.StatusOK,
		Body:       &blockingBody{},
		Request:    req,
	}, nil
}

// TestTransportPoolBalancesBeyondShardCount 验证并发超过分片数时，按在途最少策略均匀分摊，
// 且全部响应关闭后在途计数归零。
func TestTransportPoolBalancesBeyondShardCount(t *testing.T) {
	const (
		shardCount   = 4
		totalWorkers = 8
		wantPerShard = 2
		waitTimeout  = 2 * time.Second
	)
	gate := make(chan struct{})
	closeBodies := make(chan struct{})
	t.Cleanup(func() {
		select {
		case <-gate:
		default:
			close(gate)
		}
		select {
		case <-closeBodies:
		default:
			close(closeBodies)
		}
	})

	rts := make([]*barrierRoundTripper, 0, shardCount)
	pool := newTransportPool(shardCount, func() http.RoundTripper {
		rt := &barrierRoundTripper{gate: gate}
		rts = append(rts, rt)
		return rt
	})

	var returned atomic.Int32
	allReturned := make(chan struct{})
	errs := make([]error, totalWorkers)
	var wg sync.WaitGroup
	for i := 0; i < totalWorkers; i++ {
		wg.Add(1)
		go func(worker int) {
			defer wg.Done()
			req, reqErr := http.NewRequest(http.MethodPost, "https://example.com/api/house/metrics", nil)
			if reqErr != nil {
				errs[worker] = reqErr
				return
			}
			resp, rtErr := pool.RoundTrip(req)
			if rtErr != nil {
				errs[worker] = rtErr
				return
			}
			defer func() {
				if closeErr := resp.Body.Close(); closeErr != nil && errs[worker] == nil {
					errs[worker] = closeErr
				}
			}()
			if int(returned.Add(1)) == totalWorkers {
				close(allReturned)
			}
			<-closeBodies
		}(i)
	}

	deadline := time.Now().Add(waitTimeout)
	for {
		var totalInFlight int64
		for _, shard := range pool.shards {
			totalInFlight += shard.inFlight.Load()
		}
		if totalInFlight == totalWorkers {
			break
		}
		if time.Now().After(deadline) {
			close(gate)
			close(closeBodies)
			wg.Wait()
			t.Fatalf("timed out waiting for %d in-flight acquisitions, got %d", totalWorkers, totalInFlight)
		}
		time.Sleep(time.Millisecond)
	}

	for i, shard := range pool.shards {
		if got := shard.inFlight.Load(); got != wantPerShard {
			close(gate)
			close(closeBodies)
			wg.Wait()
			t.Fatalf("shard %d inFlight = %d, want %d", i, got, wantPerShard)
		}
	}

	close(gate)
	select {
	case <-allReturned:
	case <-time.After(waitTimeout):
		close(closeBodies)
		wg.Wait()
		t.Fatalf("timed out waiting for responses to return (returned=%d)", returned.Load())
	}

	// Body 尚未关闭前，在途计数仍应保持占有，防止 Body 流式读取未完成时提前抢占分片。
	for i, shard := range pool.shards {
		if got := shard.inFlight.Load(); got != wantPerShard {
			close(closeBodies)
			wg.Wait()
			t.Fatalf("shard %d inFlight before body close = %d, want %d", i, got, wantPerShard)
		}
	}

	close(closeBodies)
	wg.Wait()

	for i, workerErr := range errs {
		if workerErr != nil {
			t.Fatalf("worker %d failed: %v", i, workerErr)
		}
	}
	for i, shard := range pool.shards {
		if got := shard.inFlight.Load(); got != 0 {
			t.Fatalf("shard %d inFlight after close = %d, want 0", i, got)
		}
		if got := rts[i].maxConcurrent.Load(); got != wantPerShard {
			t.Fatalf("shard %d maxConcurrent = %d, want %d", i, got, wantPerShard)
		}
	}
}

type failingRoundTripper struct {
	err error
}

func (f failingRoundTripper) RoundTrip(*http.Request) (*http.Response, error) {
	return nil, f.err
}

// TestTransportPoolReleasesOnError 验证底层 RoundTrip 出错时立即释放在途计数，防止分片永久占满。
func TestTransportPoolReleasesOnError(t *testing.T) {
	const shardCount = 2
	boom := errors.New("dial failed")
	pool := newTransportPool(shardCount, func() http.RoundTripper {
		return failingRoundTripper{err: boom}
	})

	req, err := http.NewRequest(http.MethodGet, "https://example.com", nil)
	if err != nil {
		t.Fatalf("new request: %v", err)
	}
	resp, rtErr := pool.RoundTrip(req)
	if resp != nil && resp.Body != nil {
		_ = resp.Body.Close()
	}
	if !errors.Is(rtErr, boom) {
		t.Fatalf("RoundTrip error = %v, want %v", rtErr, boom)
	}
	for i, shard := range pool.shards {
		if got := shard.inFlight.Load(); got != 0 {
			t.Fatalf("shard %d inFlight after error = %d, want 0", i, got)
		}
	}
}

// TestTransportPoolCloseIdleConnections 验证 CloseIdleConnections 透传到全部分片。
func TestTransportPoolCloseIdleConnections(t *testing.T) {
	const shardCount = 4
	trackers := make([]*closeTracker, 0, shardCount)
	pool := newTransportPool(shardCount, func() http.RoundTripper {
		tracker := &closeTracker{}
		trackers = append(trackers, tracker)
		return tracker
	})

	pool.CloseIdleConnections()
	for i, tracker := range trackers {
		if got := tracker.closed.Load(); got != 1 {
			t.Fatalf("tracker %d CloseIdleConnections calls = %d, want 1", i, got)
		}
	}
}

// TestCloneShardTransportBudgets 验证默认分片的连接池预算、不修改默认 Transport 与 HTTP/2 保留状态。
func TestCloneShardTransportBudgets(t *testing.T) {
	base, ok := http.DefaultTransport.(*http.Transport)
	beforeIdle, beforePerHost := 0, 0
	if ok {
		beforeIdle, beforePerHost = base.MaxIdleConns, base.MaxIdleConnsPerHost
	}

	transport := cloneShardTransport(defaultTransportShards)
	wantIdle := defaultMaxIdleConns / defaultTransportShards
	wantPerHost := defaultMaxIdleConnsPerHost / defaultTransportShards
	if got := transport.MaxIdleConns; got != wantIdle {
		t.Fatalf("MaxIdleConns = %d, want %d", got, wantIdle)
	}
	if got := transport.MaxIdleConnsPerHost; got != wantPerHost {
		t.Fatalf("MaxIdleConnsPerHost = %d, want %d", got, wantPerHost)
	}
	if !transport.ForceAttemptHTTP2 {
		t.Fatal("ForceAttemptHTTP2 = false, want true so shards negotiate HTTP/2")
	}
	if ok && (base.MaxIdleConns != beforeIdle || base.MaxIdleConnsPerHost != beforePerHost) {
		t.Fatal("cloneShardTransport must not mutate http.DefaultTransport")
	}
}

// TestCloneShardTransportKeepsGlobalIdleBudgetAt16Shards 验证 16 分片时全局空闲预算仍不超过默认值：
// 每分片每主机至少保留一条空闲连接，避免分片拆细后出现 0 配额。
func TestCloneShardTransportKeepsGlobalIdleBudgetAt16Shards(t *testing.T) {
	const shards = 16
	transport := cloneShardTransport(shards)
	wantIdle := defaultMaxIdleConns / shards
	wantPerHost := defaultMaxIdleConnsPerHost / shards
	if wantPerHost < minShardIdleConnsPerHost {
		wantPerHost = minShardIdleConnsPerHost
	}
	if got := transport.MaxIdleConns; got != wantIdle {
		t.Fatalf("MaxIdleConns = %d, want %d", got, wantIdle)
	}
	if got := transport.MaxIdleConnsPerHost; got != wantPerHost {
		t.Fatalf("MaxIdleConnsPerHost = %d, want %d", got, wantPerHost)
	}
	if got := transport.MaxIdleConns * shards; got > defaultMaxIdleConns {
		t.Fatalf("global MaxIdleConns budget = %d, want <= %d", got, defaultMaxIdleConns)
	}
	if got := transport.MaxIdleConnsPerHost * shards; got > defaultMaxIdleConnsPerHost {
		t.Fatalf("global MaxIdleConnsPerHost budget = %d, want <= %d", got, defaultMaxIdleConnsPerHost)
	}
}

// TestNewFactoryUsesDefaultTransportPool 确认默认 Factory 经 retryTransport 接入 8 分片 transportPool，
// 且显式传入 BaseTransport 时保持直接透传不额外套池。
func TestNewFactoryUsesDefaultTransportPool(t *testing.T) {
	factory := NewFactory(Options{})
	t.Cleanup(func() { _ = factory.Close() })

	retryRT, ok := factory.transport.(*retryTransport)
	if !ok {
		t.Fatalf("factory.transport = %T, want *retryTransport", factory.transport)
	}
	pool, isPool := retryRT.base.(*transportPool)
	if !isPool {
		t.Fatalf("retryTransport.base = %T, want *transportPool", retryRT.base)
	}
	if got := len(pool.shards); got != defaultTransportShards {
		t.Fatalf("pool shards = %d, want %d", got, defaultTransportShards)
	}

	tracker := &closeTracker{}
	customFactory := NewFactory(Options{BaseTransport: tracker})
	customRetry, customOK := customFactory.transport.(*retryTransport)
	if !customOK {
		t.Fatalf("customFactory.transport = %T, want *retryTransport", customFactory.transport)
	}
	if customRetry.base != tracker {
		t.Fatalf("custom base = %v, want injected BaseTransport", customRetry.base)
	}
}
