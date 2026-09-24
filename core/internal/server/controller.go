// Package server 负责管理本地 gRPC 服务的宿主环境与生命周期。
//
// 职责：
// - 管理服务级生命周期状态（NotReady、Ready、Stopping、Closed）；
// - 统一仲裁服务退出流程，支持优雅退出与超时兜底强制关闭。
package server

import (
	"sync"
	"time"

	"google.golang.org/grpc"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// Controller 负责协调服务的统一关闭路径，并管理服务生命周期。
//
// 实现方式：
// 汇集 Teardown RPC、系统信号、父进程退出及 Serve 异常等所有退出来源；
// 通过 sync.Once 保证幂等执行：先触发 GracefulStop 尝试优雅退出，
// 超时后强制调用 Stop 兜底，确保进程可靠结束。
type Controller struct {
	server *grpc.Server
	grace  time.Duration
	once   sync.Once
	done   chan struct{}
	lc     *Lifecycle
}

// NewController 创建服务控制器，grace 为优雅退出的等待上限。
func NewController(g *grpc.Server, grace time.Duration) *Controller {
	return &Controller{
		server: g,
		grace:  grace,
		done:   make(chan struct{}),
		lc:     newLifecycle(),
	}
}

// Lifecycle 返回由本控制器管理的服务生命周期状态机。
func (c *Controller) Lifecycle() *Lifecycle {
	return c.lc
}

// Shutdown 触发服务关闭流程，操作具备幂等性；cause 仅用于日志标识退出来源。
func (c *Controller) Shutdown(cause string) {
	c.once.Do(func() {
		c.lc.BeginStopping()
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
			c.lc.Close()
		}()
	})
}

// Done 在关闭序列完成后关闭，供调用方等待收尾结束。
func (c *Controller) Done() <-chan struct{} {
	return c.done
}
