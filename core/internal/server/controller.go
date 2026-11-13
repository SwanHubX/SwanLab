// Package server 负责管理本地 gRPC 服务的宿主环境与生命周期。
//
// 职责：
// - 管理服务级生命周期状态（NotReady、Ready、Stopping、Closed）；
// - 统一仲裁服务退出流程，支持优雅退出与超时兜底强制关闭。
package server

import (
	"context"
	"sync"
	"time"

	"google.golang.org/grpc"

	"github.com/swanhubx/swanlab/core/internal/pkg/console"
)

// Controller 负责协调服务的统一关闭路径，并管理服务生命周期。
//
// 实现方式：
// 汇集 Teardown RPC、系统信号、父进程退出及 Serve 异常等所有退出来源；
// 通过 sync.Once 保证幂等执行：先执行 onShutdown 钩子回收活跃会话，再触发
// GracefulStop；两者共享同一个 grace deadline，超时后强制调用 Stop 兜底。
type Controller struct {
	server *grpc.Server
	grace  time.Duration
	once   sync.Once
	done   chan struct{}
	lc     *Lifecycle

	hookMu     sync.RWMutex
	onShutdown func(context.Context) error
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

// SetOnShutdown 注册进程退出钩子；钩子与 GracefulStop 共享 grace 预算。
// 必须在 Shutdown 之前调用；钩子应尊重 ctx，超时后控制器仍会进入 Stop 兜底。
func (c *Controller) SetOnShutdown(hook func(context.Context) error) {
	c.hookMu.Lock()
	c.onShutdown = hook
	c.hookMu.Unlock()
}

func (c *Controller) shutdownHook() func(context.Context) error {
	c.hookMu.RLock()
	defer c.hookMu.RUnlock()
	return c.onShutdown
}

// Shutdown 触发服务关闭流程，操作具备幂等性；cause 仅用于日志标识退出来源。
func (c *Controller) Shutdown(cause string) {
	c.once.Do(func() {
		c.lc.BeginStopping()
		console.Infof("core service shutting down (%s)", cause)
		go func() {
			defer close(c.done)

			ctx, cancel := context.WithTimeout(context.Background(), c.grace)
			defer cancel()

			if hook := c.shutdownHook(); hook != nil {
				hookDone := make(chan struct{})
				go func() {
					defer close(hookDone)
					if err := hook(ctx); err != nil {
						console.Warningf("shutdown hook failed: %v", err)
					}
				}()
				select {
				case <-hookDone:
				case <-ctx.Done():
					console.Warningf("shutdown hook exceeded grace %s", c.grace)
				}
			}

			graceful := make(chan struct{})
			go func() {
				c.server.GracefulStop()
				close(graceful)
			}()
			select {
			case <-graceful:
			case <-ctx.Done():
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
