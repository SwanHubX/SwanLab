// Package cron 提供可复用的周期任务运行时：每个任务独立协程、按间隔触发、合并积压 tick，
// 并支持注入时钟与幂等停止。
package cron

import (
	"context"
	"sync"
	"time"
)

// Clock 提供周期任务节奏；测试可注入手动推进的假时钟。
// Ticker 返回的停止函数必须由调度器执行，测试实现可以返回 no-op。
type Clock interface {
	Ticker(d time.Duration) (<-chan time.Time, func())
}

type realClock struct{}

func (realClock) Ticker(d time.Duration) (<-chan time.Time, func()) {
	t := time.NewTicker(d)
	return t.C, t.Stop
}

// Job 是一类周期任务；同一任务的执行不会重叠。
type Job struct {
	// Interval 是触发间隔，必须为正。
	Interval time.Duration
	// Run 在任务协程中同步执行；ctx 在 Stop 时取消。
	Run func(ctx context.Context)
}

// Scheduler 持有若干周期任务；Add 必须在 Start 之前，Stop 幂等且等待在途任务。
type Scheduler struct {
	clock Clock
	jobs  []Job

	mu      sync.Mutex
	started bool
	cancel  context.CancelFunc
	wg      sync.WaitGroup
}

// New 创建调度器；clock 为 nil 时使用真实定时器。
func New(clock Clock) *Scheduler {
	if clock == nil {
		clock = realClock{}
	}
	return &Scheduler{clock: clock}
}

// Add 注册任务；Interval <= 0 或 Run 为 nil 属于调用方错误，直接 panic。
func (s *Scheduler) Add(job Job) {
	if job.Interval <= 0 || job.Run == nil {
		panic("cron: Add requires positive interval and Run")
	}
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.started {
		panic("cron: Add after Start")
	}
	s.jobs = append(s.jobs, job)
}

// Start 为每个任务启动一个协程；重复调用无效。
func (s *Scheduler) Start(parent context.Context) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.started {
		return
	}
	s.started = true
	if parent == nil {
		parent = context.Background()
	}
	ctx, cancel := context.WithCancel(parent)
	s.cancel = cancel
	for _, job := range s.jobs {
		s.wg.Add(1)
		go func(job Job) {
			defer s.wg.Done()
			s.loop(ctx, job)
		}(job)
	}
}

func (s *Scheduler) loop(ctx context.Context, job Job) {
	ticks, stop := s.clock.Ticker(job.Interval)
	defer stop()
	for {
		select {
		case <-ctx.Done():
			return
		case _, ok := <-ticks:
			if !ok {
				return
			}
		}
		// 合并等待期间积压的 tick：慢任务不追帧，只按最新一次触发执行。
	drain:
		for {
			select {
			case <-ticks:
			default:
				break drain
			}
		}
		job.Run(ctx)
	}
}

// Stop 取消全部任务并等待在途执行结束；可重复调用。
func (s *Scheduler) Stop() {
	s.mu.Lock()
	if s.cancel != nil {
		s.cancel()
	}
	s.mu.Unlock()
	s.wg.Wait()
}
