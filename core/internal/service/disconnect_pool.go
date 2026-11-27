package service

import (
	"context"
	"sync"
	"time"
)

const defaultDisconnectGrace = 30 * time.Second

// disconnectPool 管理传输断连后的异步会话清理协程，并在服务收尾时阻塞等待或取消。
type disconnectPool struct {
	mu     sync.Mutex
	closed bool
	grace  time.Duration
	ctx    context.Context
	cancel context.CancelFunc
	wg     sync.WaitGroup
}

func newDisconnectPool(grace time.Duration) *disconnectPool {
	if grace <= 0 {
		grace = defaultDisconnectGrace
	}
	ctx, cancel := context.WithCancel(context.Background())
	return &disconnectPool{
		grace:  grace,
		ctx:    ctx,
		cancel: cancel,
	}
}

func (p *disconnectPool) isClosed() bool {
	p.mu.Lock()
	defer p.mu.Unlock()
	return p.closed
}

func (p *disconnectPool) goRun(fn func(context.Context)) {
	p.mu.Lock()
	if p.closed {
		p.mu.Unlock()
		return
	}
	p.wg.Add(1)
	grace := p.grace
	p.mu.Unlock()

	go func() {
		defer p.wg.Done()
		ctx, cancel := context.WithTimeout(p.ctx, grace)
		defer cancel()
		fn(ctx)
	}()
}

func (p *disconnectPool) wait(ctx context.Context) error {
	p.mu.Lock()
	p.closed = true
	p.mu.Unlock()

	done := make(chan struct{})
	go func() {
		p.wg.Wait()
		close(done)
	}()

	select {
	case <-done:
		p.cancel()
		return nil
	case <-ctx.Done():
		p.cancel()
		<-done
		return ctx.Err()
	}
}
