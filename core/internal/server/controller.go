package server

import (
	"sync"
	"time"

	"google.golang.org/grpc"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// Controller 仲裁服务的统一关闭路径。
//
// Teardown RPC、SIGINT/SIGTERM、父进程退出通知与 Serve 异常都汇入同一条
// 收尾序列：先 GracefulStop 等待在途请求完成，超过 grace 时限后强制 Stop，
// 并保证 Serve 一定返回。Shutdown 幂等，多次触发只执行一次。
type Controller struct {
	server *grpc.Server
	grace  time.Duration
	once   sync.Once
	done   chan struct{}
}

// NewController 包装一个 gRPC Server，grace 为优雅关闭的等待上限。
func NewController(g *grpc.Server, grace time.Duration) *Controller {
	return &Controller{
		server: g,
		grace:  grace,
		done:   make(chan struct{}),
	}
}

// Shutdown 幂等触发关闭；cause 仅用于日志，标识关闭来源。
func (c *Controller) Shutdown(cause string) {
	c.once.Do(func() {
		console.Infof("core service shutting down (%s)", cause)
		go func() {
			defer close(c.done)
			graceful := make(chan struct{})
			go func() {
				c.server.GracefulStop()
				close(graceful)
			}()
			select {
			case <-graceful:
			case <-time.After(c.grace):
				c.server.Stop()
				<-graceful
			}
		}()
	})
}

// Done 在关闭序列完成后关闭，供调用方等待收尾结束。
func (c *Controller) Done() <-chan struct{} {
	return c.done
}
