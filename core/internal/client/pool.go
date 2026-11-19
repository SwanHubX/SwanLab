package client

import (
	"io"
	"net/http"
	"sync"
	"sync/atomic"
)

const (
	// defaultTransportShards 是共享客户端的默认连接分片数。
	// Go http.Transport 在 HTTP/2 下会将同 host 并发请求全部复用到单条 TCP 连接，
	// 拆分为多个独立分片可在保留 HTTP/2 协议的前提下展开多条并行连接。
	defaultTransportShards = 8

	// 分片拆小时仍为每片保留最小空闲池额度，避免高频并发下频繁关闭重建连接。
	defaultMinShardIdleConns        = 4
	defaultMinShardIdleConnsPerHost = 2
	minPoolShards                   = 1
)

type transportShard struct {
	rt       http.RoundTripper
	inFlight atomic.Int64
}

func (s *transportShard) release() {
	s.inFlight.Add(-1)
}

// transportPool 把并发请求分散到多个独立的 RoundTripper 分片，解除 HTTP/2 单连接上限。
// 低并发下优先抢占最靠前的空闲分片（让单会话串行流量持续复用分片 0 的已建连接），
// 仅当分片忙时才向外扩张或按在途最少策略均衡分摊。
type transportPool struct {
	shards []*transportShard
	cursor atomic.Uint64
}

var _ http.RoundTripper = (*transportPool)(nil)

type releaseBody struct {
	io.ReadCloser
	once    sync.Once
	release func()
}

func (b *releaseBody) Close() error {
	defer b.once.Do(b.release)
	return b.ReadCloser.Close()
}

// cloneShardTransport 构造单个分片的 HTTP Transport，保留 HTTP/2 协商能力，
// 并按分片数均分空闲连接预算，使全局总配额仍维持在默认上限附近。
func cloneShardTransport(shards int) *http.Transport {
	transport := cloneDefaultTransport()
	if shards <= minPoolShards {
		return transport
	}
	transport.MaxIdleConns = max(defaultMinShardIdleConns, defaultMaxIdleConns/shards)
	transport.MaxIdleConnsPerHost = max(defaultMinShardIdleConnsPerHost, defaultMaxIdleConnsPerHost/shards)
	return transport
}

func newDefaultTransportPool(shards int) *transportPool {
	return newTransportPool(shards, func() http.RoundTripper {
		return cloneShardTransport(shards)
	})
}

func newTransportPool(shards int, makeShard func() http.RoundTripper) *transportPool {
	if shards < minPoolShards {
		shards = minPoolShards
	}
	if makeShard == nil {
		makeShard = func() http.RoundTripper {
			return cloneShardTransport(shards)
		}
	}
	items := make([]*transportShard, shards)
	for i := range items {
		items[i] = &transportShard{rt: makeShard()}
	}
	return &transportPool{shards: items}
}

func (p *transportPool) acquire() *transportShard {
	for _, shard := range p.shards {
		if shard.inFlight.CompareAndSwap(0, 1) {
			return shard
		}
	}

	n := len(p.shards)
	start := int((p.cursor.Add(1) - 1) % uint64(n))
	best := p.shards[start]
	bestLoad := best.inFlight.Load()
	for offset := 1; offset < n; offset++ {
		candidate := p.shards[(start+offset)%n]
		if load := candidate.inFlight.Load(); load < bestLoad {
			best = candidate
			bestLoad = load
		}
	}
	best.inFlight.Add(1)
	return best
}

// RoundTrip 选择一个分片发送请求，并在响应体关闭（或请求失败）时归还在途计数。
func (p *transportPool) RoundTrip(req *http.Request) (*http.Response, error) {
	shard := p.acquire()
	releasedByBody := false
	defer func() {
		if !releasedByBody {
			shard.release()
		}
	}()

	resp, err := shard.rt.RoundTrip(req)
	if err != nil || resp == nil || resp.Body == nil {
		return resp, err
	}
	resp.Body = &releaseBody{
		ReadCloser: resp.Body,
		release:    shard.release,
	}
	releasedByBody = true
	return resp, nil
}

// CloseIdleConnections 释放全部分片的空闲连接。
func (p *transportPool) CloseIdleConnections() {
	for _, shard := range p.shards {
		if closer, ok := shard.rt.(interface{ CloseIdleConnections() }); ok {
			closer.CloseIdleConnections()
		}
	}
}
