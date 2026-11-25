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
	"github.com/swanhubx/swanlab/core/internal/manager/backfill"
	"github.com/swanhubx/swanlab/core/internal/manager/checkpoint"
	"github.com/swanhubx/swanlab/core/internal/manager/heartbeat"
	"github.com/swanhubx/swanlab/core/internal/manager/progress"
	"github.com/swanhubx/swanlab/core/internal/manager/recovery"
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
	// queueFactor 是通道容量相对单批记录数的倍数，提供更高的突发写入缓冲余量。
	queueFactor = 8
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
	// finishMu 串行化收尾序列：DeliverRunFinish、后台在线排空、Confirm 续跑与 shutdown 共用。
	finishMu sync.Mutex
	// finishCommitted 标记 FinishRecord 已完成提交（online 在本地存储降级时跳过），避免重试时重复写入。
	finishCommitted bool
	// localFinished 标记本地收尾已完成，供重复 DeliverRunFinish 快速返回而不被后台补投阻塞。
	localFinished atomic.Bool
	// uploadDrainPrepared 标记在线补投、崩溃日志入队与 Intake.Close 已完成。
	uploadDrainPrepared bool
	// pendingLog 是崩溃 error_log：由 mu 保护，只上传不落盘，暂时入队失败时保留供重试。
	pendingLog *recordv1.Record

	// shutdownCtx 在会话退出或关闭时取消，用于立即打断后台上传和结束重试。
	shutdownCtx    context.Context
	shutdownCancel context.CancelFunc
	drainMu        sync.Mutex
	drainStarted   bool
	drainDone      chan struct{}

	state atomic.Int32

	// drained 在会话完成上传时关闭；failed 在结束流程失败时关闭。
	drained chan struct{}
	failed  chan struct{}

	drainFn func(context.Context) error
	cloud   atomic.Pointer[cloudResources]

	// 云端连接状态和后台任务的生命周期，由 cloudMu 保护。
	cloudMu      sync.Mutex
	cloudState   operationv1.CloudState
	cloudPath    string
	cloudErr     error
	attachDone   chan struct{}
	attachOnce   sync.Once
	attachCancel context.CancelFunc
	attachWG     sync.WaitGroup

	finishErr error

	// 装配期依赖；Start 失败回滚后可能只完成部分初始化，使用前必须判空。
	ctx     context.Context
	cancel  context.CancelFunc
	cfg     config.Settings
	run     *runctx.RunContext
	store   *store.AsyncWriter
	intake  *buffer.Intake
	budget  *budget.Budget
	tracker *uploader.Tracker
	watcher watcher.Watcher

	// 指标与编号状态由 Session.mu 保护；Resume 恢复的 epoch 作为日志起始行号。
	metrics *metrics.RunMetrics
	counter *counter.Counter
	epoch   *counter.Counter
	// finishRecord 保存最近的结束记录，Confirm 上报 stop_experiment 时读取终态。
	finishRecord *runv1.FinishRecord

	storeDegraded atomic.Bool

	// recovery 控制 per-session 的本地失效补传协程。
	recovery *recovery.Pump

	// 补投：backfill 管理暂停区间的后台重放与停机补完。
	backfill *backfill.Worker

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

	progress   *progress.Tracker
	checkpoint *checkpoint.Writer
}

// cloudResources 保存上传器、心跳和停止实验的函数；创建后只读，通过原子指针保存。
type cloudResources struct {
	uploader  *uploader.Uploader
	heartbeat *heartbeat.Heartbeat
	releaseFn func(context.Context) error
}

func (s *Session) setCloudResources(res *cloudResources) {
	s.cloud.Store(res)
}

func (s *Session) uploader() *uploader.Uploader {
	if res := s.cloud.Load(); res != nil {
		return res.uploader
	}
	return nil
}

func (s *Session) heartbeat() *heartbeat.Heartbeat {
	if res := s.cloud.Load(); res != nil {
		return res.heartbeat
	}
	return nil
}

func (s *Session) releaseFn() func(context.Context) error {
	if res := s.cloud.Load(); res != nil {
		return res.releaseFn
	}
	return nil
}

func (s *Session) isCloudAttached() bool {
	return s.uploader() != nil
}

func (s *Session) initDeferredAttach() {
	s.cloudMu.Lock()
	defer s.cloudMu.Unlock()
	s.cloudState = operationv1.CloudState_CLOUD_STATE_PENDING
	s.cloudPath = ""
	s.cloudErr = nil
	s.attachDone = make(chan struct{})
}

func (s *Session) cloudStatus() (operationv1.CloudState, string, error) {
	s.cloudMu.Lock()
	defer s.cloudMu.Unlock()
	return s.cloudState, s.cloudPath, s.cloudErr
}

func (s *Session) closeAttachDone() {
	s.cloudMu.Lock()
	done := s.attachDone
	s.cloudMu.Unlock()
	if done != nil {
		s.attachOnce.Do(func() { close(done) })
	}
}

func (s *Session) markCloudAttached(path string) {
	s.cloudMu.Lock()
	s.cloudState = operationv1.CloudState_CLOUD_STATE_ATTACHED
	s.cloudPath = path
	s.cloudErr = nil
	s.cloudMu.Unlock()
	s.closeAttachDone()
}

// markAttachStopped 在连接任务因会话结束退出时关闭 attachDone。
// 不修改 cloudState：既不能标记为成功，也不能误报为永久失败。
func (s *Session) markAttachStopped() {
	s.closeAttachDone()
}

func (s *Session) isPermanentAttachFailure() bool {
	state, _, _ := s.cloudStatus()
	return state == operationv1.CloudState_CLOUD_STATE_FAILED
}

// localCopyFailure 按照本地副本状态拼接失败文案；suffix 只用于本地副本不完整的两种情况。
func (s *Session) localCopyFailure(prefix, suffix string) string {
	switch {
	case s.cfg.SkipStore:
		return fmt.Sprintf("%s, and no local copy was kept (skip_store)%s.", prefix, suffix)
	case s.storeDegraded.Load():
		return fmt.Sprintf("%s, and the local copy is incomplete%s.", prefix, suffix)
	default:
		return fmt.Sprintf("%s, but it has been saved locally.", prefix)
	}
}

func (s *Session) attachFailureConfirmMessage(attachErr error) string {
	reason := "connection failed"
	if attachErr != nil {
		reason = attachErr.Error()
	}
	return s.localCopyFailure(fmt.Sprintf("Failed to create run on SwanLab online service (%s)", reason), "")
}

// noteAttachError 记录最近一次可重试的连接错误，供进度界面显示等待原因。
// 状态保持 PENDING，连接成功后由 markCloudAttached 清空。
func (s *Session) noteAttachError(err error) {
	if err == nil {
		return
	}
	s.cloudMu.Lock()
	if s.cloudState == operationv1.CloudState_CLOUD_STATE_PENDING {
		s.cloudErr = err
	}
	s.cloudMu.Unlock()
}

func (s *Session) markCloudFailed(err error) {
	s.cloudMu.Lock()
	s.cloudState = operationv1.CloudState_CLOUD_STATE_FAILED
	s.cloudErr = err
	s.cloudMu.Unlock()
	s.closeAttachDone()
}

func newSession(opts sessionOptions) *Session {
	shutdownCtx, shutdownCancel := context.WithCancel(context.Background())
	s := &Session{
		drained:        make(chan struct{}),
		failed:         make(chan struct{}),
		drainFn:        opts.drainFn,
		healthWake:     make(chan struct{}, 1),
		shutdownCtx:    shutdownCtx,
		shutdownCancel: shutdownCancel,
	}
	if opts.releaseFn != nil {
		s.setCloudResources(&cloudResources{releaseFn: opts.releaseFn})
	}
	s.state.Store(int32(RunStateRunning))
	return s
}

// withShutdownCancel 创建一个同时受父 ctx 和会话关闭信号控制的子 ctx。
func (s *Session) withShutdownCancel(parent context.Context) (context.Context, context.CancelFunc) {
	if parent == nil {
		parent = context.Background()
	}
	ctx, cancel := context.WithCancel(parent)
	if s.shutdownCtx == nil {
		return ctx, cancel
	}
	stop := context.AfterFunc(s.shutdownCtx, cancel)
	return ctx, func() {
		stop()
		cancel()
	}
}

// startOnlineDrain 在本地结束后只启动一次后台上传任务。
func (s *Session) startOnlineDrain() {
	s.drainMu.Lock()
	if s.drainStarted || (s.shutdownCtx != nil && s.shutdownCtx.Err() != nil) {
		s.drainMu.Unlock()
		return
	}
	s.drainStarted = true
	done := make(chan struct{})
	s.drainDone = done
	s.drainMu.Unlock()

	go func() {
		defer close(done)
		drainCtx, cancel := s.withShutdownCancel(s.ctx)
		defer cancel()
		s.runOnlineDrain(drainCtx)
	}()
}

// stopOnlineDrain 通知会话关闭并等待后台上传任务退出；调用时不能持有 finishMu。
func (s *Session) stopOnlineDrain() {
	s.drainMu.Lock()
	if s.shutdownCancel != nil {
		s.shutdownCancel()
	}
	done := s.drainDone
	s.drainMu.Unlock()
	if done != nil {
		<-done
	}
}

func (s *Session) crashLog() *recordv1.Record {
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.pendingLog
}

func (s *Session) clearCrashLog() {
	s.mu.Lock()
	s.pendingLog = nil
	s.mu.Unlock()
}

func (s *Session) hasPendingLog() bool {
	return s.crashLog() != nil
}

// beginDrain 把会话从 RUNNING 切到 DRAINING；重复调用无副作用。
// 已 FAILED 时返回原始失败，调用方不能当作成功。
func (s *Session) beginDrain() error {
	_, err := s.beginDrainOwned()
	return err
}

// beginDrainOwned 切换状态，并报告本次调用是否真正完成了 RUNNING → DRAINING。
// 只有完成切换的调用方才能执行第一次结束流程；FAILED 的后续处理由 finish/resumeFinish 负责。
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

// markDrained 由上传完成回调调用；只允许 DRAINING → DRAINED。
func (s *Session) markDrained() bool {
	if !s.state.CompareAndSwap(int32(RunStateDraining), int32(RunStateDrained)) {
		return false
	}
	s.setTrackerFinished()
	close(s.drained)
	return true
}

// markRecovered 在 Confirm 重试成功后把 FAILED 恢复为 DRAINED。
func (s *Session) markRecovered() bool {
	if !s.state.CompareAndSwap(int32(RunStateFailed), int32(RunStateDrained)) {
		return false
	}
	s.setTrackerFinished()
	close(s.drained)
	return true
}

// setTrackerFinished 让进度轮询在终态停止等待；FAILED 也要设置。
func (s *Session) setTrackerFinished() {
	if s.tracker != nil {
		s.tracker.SetState(operationv1.CoreState_CORE_STATE_FINISHED)
	}
}

// fail 记录结束失败并把会话置为 FAILED；DRAINING 和 FAILED 都允许更新原因，
// 便于 Confirm 再次失败时返回最新信息。
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

// failure 返回结束失败信息。
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
	// 保留原始原因，同时让调用方可以用 ErrRunFailed 判断“结束失败”。
	return fmt.Errorf("%w: %w", ErrRunFailed, s.finishErr)
}

// stateSnapshot 返回当前会话状态的并发安全快照。
func (s *Session) stateSnapshot() RunState {
	return RunState(s.state.Load())
}

// waitDrained 在锁外等待会话完成上传，返回 nil 表示可以进入释放阶段。
// FAILED 也返回 nil：Confirm 会在 confirmMu 内重新执行被中断的结束步骤。
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

	// 结束失败后 drainFn 不会返回（intake 未关闭，上传器无法完成），
	// 因此在独立协程中等待，让 FAILED 能提前唤醒调用方；提前返回时取消子 ctx。
	waitCtx, cancel := context.WithCancel(ctx)
	defer cancel()
	drained := make(chan error, 1)
	go func() { drained <- s.drainFn(waitCtx) }()

	select {
	case err := <-drained:
		// Wait 返回 nil 时 DRAINED 已设置：OnDrained 先于 Done 关闭执行。
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
		// 与上面一致：把重新执行的机会留给 Confirm。
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
	cloudState, cloudPath, cloudErr := s.cloudStatus()
	if cloudState != operationv1.CloudState_CLOUD_STATE_UNSPECIFIED {
		errText := ""
		if cloudErr != nil {
			errText = cloudErr.Error()
		}
		stats.Cloud = &operationv1.CloudStatus{
			State: cloudState,
			Path:  cloudPath,
			Error: errText,
		}
	}
	up := s.uploader()
	if up == nil && s.intake == nil && s.budget == nil {
		return stats
	}
	upload := &operationv1.UploadStatus{}
	if up != nil {
		upload.PendingRecords = int64(up.Stats().Pending)
		upload.Paused = up.Suspended()
		upload.Retrying = up.InFailure()
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
