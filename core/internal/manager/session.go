// Package manager 负责维护 swanlab-core 的 Run 会话与业务编排逻辑。
//
// 该包不依赖 gRPC 服务契约，但直接消费数据面 proto（run/data/column/save/terminal），
// 专职管理 Run 级状态机、会话装配与云端编排。
package manager

import (
	"context"
	"errors"
	"fmt"
	"sync"
	"sync/atomic"
	"time"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/cron/health"
	"github.com/swanhubx/swanlab/core/internal/manager/heartbeat"
	"github.com/swanhubx/swanlab/core/internal/manager/progress"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

const (
	// rollbackTimeout 限制启动失败回滚时等待 intake 排空的时间。
	rollbackTimeout = time.Second
	// queueFactor 是通道容量相对单批记录数的倍数：自动列最多让一批记录翻倍。
	queueFactor = 2
	// finishUnknownError 与 Python 对齐：FinishRecord 未带 error 时的默认错误文案。
	finishUnknownError = "run failed with unknown error while finish_record.error is not set"
)

// errFinishRecord 标记 FinishRecord 未写入本地副本：Confirm 必须重试提交，
// 不能在缺少终态记录的情况下上报成功。
var errFinishRecord = errors.New("finish record was not persisted")

// RunState 表示单个实验会话的内部生命周期状态。
type RunState int32

const (
	// RunStateRunning 会话正在运行中，接受公共提交。
	RunStateRunning RunState = iota
	// RunStateDraining 会话已停止接受公共提交，正在等待上传与落盘排空。
	RunStateDraining
	// RunStateDrained 上传已排空，确认终态记录已提交后即可释放。
	RunStateDrained
	// RunStateFailed 收尾序列失败，会话不能按成功路径释放。
	RunStateFailed
)

// String 返回状态的可读名称，便于日志与测试失败信息定位。
func (s RunState) String() string {
	switch s {
	case RunStateRunning:
		return "running"
	case RunStateDraining:
		return "draining"
	case RunStateDrained:
		return "drained"
	case RunStateFailed:
		return "failed"
	default:
		return "unknown"
	}
}

// StartResult 是 DeliverRunStart 的返回值，供 service 填充 gRPC 响应。
type StartResult struct {
	SessionID        string
	Run              *runv1.StartRecord
	Path             string
	Name             string
	GlobalStep       int64
	GlobalSystemStep int64
	NewExperiment    bool
}

// sessionOptions 注入会话的排空与释放钩子。
// releaseFn 由装配期接到 stop_experiment；drainFn 接到 uploader.Wait。
type sessionOptions struct {
	drainFn   func(context.Context) error
	releaseFn func(context.Context) error
}

// sessionInit 是会话装配所需的已认证依赖。
type sessionInit struct {
	cfg          config.Settings
	run          *runctx.RunContext
	http         *client.Client
	api          *api.Client
	workspace    string
	project      string
	projectID    string
	experimentID string
	// healthReg 是进程级探针注册表；healthOpts 描述本会话所属站点的探测配置。
	healthReg  *health.Registry
	healthOpts health.SiteOptions
	// resume 非空表示基于已有 run-*.swanlab 的崩溃续传，open 以追加写和重放模式装配。
	resume *metrics.RebuildResult
}

// Session 保存单个 session_id 对应的会话状态与生命周期数据。
type Session struct {
	// mu 保护状态切换与失败信息，与公共提交共用同一把锁。
	mu sync.Mutex
	// confirmMu 只串行化「释放」阶段，不覆盖排空等待。
	confirmMu sync.Mutex
	// finishMu 串行化收尾尾部：重复的 DeliverRunFinish、Confirm 续跑与进程
	// 退出时的 shutdown 共用同一把锁，避免关闭本地通道与写入并发。
	finishMu sync.Mutex
	// finishCommitted 标记 FinishRecord 已完成提交（online 在本地存储降级时跳过），避免重试时重复写入。
	finishCommitted bool
	// pendingLog 是崩溃 error_log：只上传不落盘，暂时入队失败时保留供重试。
	pendingLog *recordv1.Record

	state atomic.Int32

	// drained 在 DRAINING → DRAINED 时关闭；failed 在进入 FAILED 时关闭。
	drained chan struct{}
	failed  chan struct{}

	drainFn   func(context.Context) error
	releaseFn func(context.Context) error

	finishErr error

	// 装配期依赖；Start 失败回滚后可能只完成部分初始化，使用前必须判空。
	ctx       context.Context
	cancel    context.CancelFunc
	cfg       config.Settings
	run       *runctx.RunContext
	store     *store.AsyncWriter
	intake    *buffer.Intake
	uploader  *uploader.Uploader
	budget    *budget.Budget
	tracker   *uploader.Tracker
	watcher   watcher.Watcher
	heartbeat *heartbeat.Heartbeat

	// 指标与编号状态由 Session.mu 保护；Resume 恢复的 epoch 作为日志起始行号。
	metrics *metrics.RunMetrics
	counter *counter.Counter
	epoch   *counter.Counter
	// finishRecord 保存最近一次收尾记录，Confirm 上报 stop_experiment 时读取终态。
	finishRecord *runv1.FinishRecord

	storeDegraded atomic.Bool

	// recoveryStop / recoveryDone 控制 per-session 的本地失效补传协程。
	recoveryStop chan struct{}
	recoveryDone chan struct{}
	recoveryOnce sync.Once

	// 补投：backfillWake 唤醒后台协程，backfillCancel 取消当前尝试。
	// backfillState 只在补投协程或停止补投后的收尾路径中访问。
	backfillWake   chan struct{}
	backfillStop   chan struct{}
	backfillDone   chan struct{}
	backfillOnce   sync.Once
	backfillMu     sync.Mutex
	backfillCancel context.CancelFunc
	backfillState  backfillState

	// 站点健康：healthReg 持有共享引用，订阅与站点由 closeHealth 幂等释放。
	healthReg         *health.Registry
	site              *health.Site
	healthUpdates     <-chan struct{}
	healthUnsubscribe func()
	healthStop        context.CancelFunc
	healthWG          sync.WaitGroup
	healthOnce        sync.Once
	healthWake        chan struct{}

	runFile        string
	runFileCreated bool

	progress           *progress.Tracker
	checkpointFile     string
	checkpointInterval time.Duration
	checkpointTicked   func()
	checkpointStop     chan struct{}
	checkpointDone     chan struct{}
	checkpointOnce     sync.Once
	checkpointMu       sync.Mutex
	checkpointWarned   bool
	checkpointLastSeq  uint64
	checkpointHasWrite bool
}

func newSession(opts sessionOptions) *Session {
	s := &Session{
		drained:    make(chan struct{}),
		failed:     make(chan struct{}),
		drainFn:    opts.drainFn,
		releaseFn:  opts.releaseFn,
		healthWake: make(chan struct{}, 1),
	}
	s.state.Store(int32(RunStateRunning))
	return s
}

// beginDrain 把会话从 RUNNING 切到 DRAINING；重复调用幂等。
// 已 FAILED 时返回原始失败，调用方不得把它当作成功。
func (s *Session) beginDrain() error {
	_, err := s.beginDrainOwned()
	return err
}

// beginDrainOwned 切换排空状态，并报告本次调用是否真正完成了 RUNNING → DRAINING。
// 只有赢得切换的调用方才能执行首次收尾；FAILED 之后的续跑由 finish/resumeFinish 处理。
func (s *Session) beginDrainOwned() (bool, error) {
	s.mu.Lock()
	defer s.mu.Unlock()
	switch RunState(s.state.Load()) {
	case RunStateRunning:
		s.state.Store(int32(RunStateDraining))
		return true, nil
	case RunStateDraining, RunStateDrained:
		return false, nil
	case RunStateFailed:
		return false, s.failureLocked()
	default:
		return false, ErrRunFailed
	}
}

// markDrained 由排空完成回调调用；只允许 DRAINING → DRAINED。
func (s *Session) markDrained() bool {
	if !s.state.CompareAndSwap(int32(RunStateDraining), int32(RunStateDrained)) {
		return false
	}
	s.setTrackerFinished()
	close(s.drained)
	return true
}

// markRecovered 在 Confirm 重试收尾成功后把 FAILED 拉回 DRAINED。
func (s *Session) markRecovered() bool {
	if !s.state.CompareAndSwap(int32(RunStateFailed), int32(RunStateDrained)) {
		return false
	}
	s.setTrackerFinished()
	close(s.drained)
	return true
}

// setTrackerFinished 让进度轮询在终态停止等待；FAILED 也要置位。
func (s *Session) setTrackerFinished() {
	if s.tracker != nil {
		s.tracker.SetState(operationv1.CoreState_CORE_STATE_FINISHED)
	}
}

// fail 记录收尾失败并把会话置为 FAILED；DRAINING 与 FAILED 都允许更新失败原因，
// 便于 Confirm 重试再次失败时返回最新的错误文案。
func (s *Session) fail(err error) {
	if err == nil {
		err = ErrRunFailed
	}
	s.mu.Lock()
	defer s.mu.Unlock()
	state := RunState(s.state.Load())
	if state != RunStateDraining && state != RunStateFailed {
		return
	}
	s.finishErr = err
	s.setTrackerFinished()
	if state == RunStateDraining {
		s.state.Store(int32(RunStateFailed))
		close(s.failed)
	}
}

// failure 返回收尾失败信息。
func (s *Session) failure() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.failureLocked()
}

func (s *Session) failureLocked() error {
	if s.finishErr == nil {
		return ErrRunFailed
	}
	if errors.Is(s.finishErr, ErrRunFailed) {
		return s.finishErr
	}
	// 保留原始原因，同时让调用方可以稳定地用 ErrRunFailed 判断“收尾失败”。
	return fmt.Errorf("%w: %w", ErrRunFailed, s.finishErr)
}

// stateSnapshot 返回当前会话状态的并发安全快照。
func (s *Session) stateSnapshot() RunState {
	return RunState(s.state.Load())
}

// waitDrained 在锁外等待会话排空，返回 nil 表示可以进入释放阶段。
// FAILED 也返回 nil：Confirm 会在 confirmMu 内重跑被中断的收尾步骤。
func (s *Session) waitDrained(ctx context.Context) error {
	switch RunState(s.state.Load()) {
	case RunStateRunning:
		return ErrRunNotFinished
	case RunStateDrained, RunStateFailed:
		return nil
	}

	if s.drainFn == nil {
		s.markDrained()
		return nil
	}

	// drainFn 在收尾失败后不会返回（intake 未关，uploader 无法排空），
	// 因此在独立协程中等待，让 FAILED 能提前唤醒调用方；提前返回时取消子 ctx。
	waitCtx, cancel := context.WithCancel(ctx)
	defer cancel()
	drained := make(chan error, 1)
	go func() { drained <- s.drainFn(waitCtx) }()

	select {
	case err := <-drained:
		// Wait 返回 nil 时 DRAINED 已置位：OnDrained 先于 Done 关闭执行。
		if err == nil {
			return nil
		}
		select {
		case <-s.failed:
			return nil
		case <-ctx.Done():
			return ctx.Err()
		default:
			return err
		}
	case <-s.drained:
		return nil
	case <-s.failed:
		// 与上面一致：把重跑机会留给 Confirm。
		return nil
	case <-ctx.Done():
		return ctx.Err()
	}
}

// snapshot 返回上传进度与投递状态；装配完成前 tracker 为空，返回空统计。
func (s *Session) snapshot() *operationv1.OperationStats {
	stats := &operationv1.OperationStats{}
	if s.tracker != nil {
		stats = s.tracker.Snapshot()
	}
	if s.uploader == nil && s.intake == nil && s.budget == nil {
		return stats
	}
	upload := &operationv1.UploadStatus{}
	if s.uploader != nil {
		upload.PendingRecords = int64(s.uploader.Stats().Pending)
		upload.Paused = s.uploader.Suspended()
	}
	if s.intake != nil {
		if paused, seq := s.intake.PauseState(); paused {
			upload.Paused = true
			upload.PauseSeq = seq
		}
	}
	if s.budget != nil {
		_, usedBytes := s.budget.Usage()
		upload.BufferedBytes = usedBytes
	}
	stats.Upload = upload
	return stats
}
