// Package manager 负责维护 swanlab-core 的 Run 会话与业务编排逻辑。
//
// 该包不依赖 gRPC 服务契约，但直接消费数据面 proto（run/data/column/save/terminal），
// 专职管理 Run 级状态机、会话装配与云端编排。
package manager

import (
	"context"
	"errors"
	"fmt"
	"os"
	"sync"
	"sync/atomic"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/pkg/watcher"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	"github.com/swanhubx/swanlab/core/internal/uploader/sender"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

const (
	// rollbackTimeout 限制启动失败回滚时等待 intake 排空的时间。
	rollbackTimeout = time.Second
	// queueFactor 是通道容量相对单批记录数的倍数：自动列最多让一批记录翻倍。
	queueFactor = 2
	// finishUnknownError 与 Python 对齐：FinishRecord 未带 error 时的默认错误文案。
	finishUnknownError = "run failed with unknown error while finish_record.error is not set"
)

// RunState 表示单个实验会话的内部生命周期状态。
type RunState int32

const (
	// RunStateRunning 会话正在运行中，接受公共提交。
	RunStateRunning RunState = iota
	// RunStateDraining 会话已停止接受公共提交，正在等待上传与落盘排空。
	RunStateDraining
	// RunStateDrained 会话已排空，可以安全释放。
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
// C6 会把 releaseFn 接到 stop_experiment；drainFn 由装配期接到 uploader.Wait。
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
}

// Session 保存单个 session_id 对应的会话状态与生命周期数据。
type Session struct {
	// mu 保护状态切换与失败信息，与公共提交共用同一把锁。
	mu sync.Mutex
	// confirmMu 只串行化「释放」阶段，不覆盖排空等待。
	confirmMu sync.Mutex

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
	tracker   *uploader.Tracker
	watcher   watcher.Watcher
	heartbeat *heartbeat

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

	runFile        string
	runFileCreated bool
}

func newSession(opts sessionOptions) *Session {
	s := &Session{
		drained:   make(chan struct{}),
		failed:    make(chan struct{}),
		drainFn:   opts.drainFn,
		releaseFn: opts.releaseFn,
	}
	s.state.Store(int32(RunStateRunning))
	return s
}

// DeliverRunStart 完成云端项目/实验编排，装配会话资源并投递 StartRecord。
//
// 任一阶段失败都会回滚已创建的资源（包括本次排他创建的 run-*.swanlab），
// 注册表只在装配全部成功后插入，因此失败不会留下不可用会话。
func (m *Manager) DeliverRunStart(
	ctx context.Context,
	cfg config.Settings,
	rec *runv1.StartRecord,
) (*StartResult, error) {
	if rec == nil {
		return nil, fmt.Errorf("%w: start record must not be nil", ErrInvalidStartRecord)
	}
	run, err := runctx.New(cfg)
	if err != nil {
		return nil, fmt.Errorf("%w: %w", ErrInvalidSettings, err)
	}
	if !m.hasCapacity() {
		return nil, fmt.Errorf("%w: limit %d", ErrTooManySessions, maxSessions)
	}
	resumeMode, err := resumeModeFromProto(rec.GetResume())
	if err != nil {
		return nil, fmt.Errorf("%w: %w", ErrInvalidStartRecord, err)
	}

	httpClient, err := m.factory.Client(ctx, cfg.APIHost, cfg.APIKey)
	if err != nil {
		return nil, fmt.Errorf("start: authenticate: %w", err)
	}
	apiClient := api.New(httpClient)

	project, err := apiClient.GetOrCreateProject(ctx, rec.GetWorkspace(), rec.GetProject(), rec.GetPublic())
	if err != nil {
		return nil, fmt.Errorf("start: get project %q: %w", rec.GetProject(), err)
	}
	// 与 Python 一致：项目查询使用 record.workspace，后续实验 owner 以项目所属组为准。
	workspace := firstNonEmpty(project.Group.Username, rec.GetWorkspace(), apiClient.Username())
	projectName := firstNonEmpty(project.Name, rec.GetProject())
	if workspace == "" || projectName == "" || project.CUID == "" || project.Path == "" {
		return nil, fmt.Errorf(
			"%w: incomplete project metadata (workspace=%q project=%q cuid=%q path=%q)",
			ErrInvalidStartRecord, workspace, projectName, project.CUID, project.Path)
	}

	name := rec.GetName()
	if name == "" {
		name = builder.GenerateName(project.Count.Experiments)
	}
	color := rec.GetColor()
	if color == "" {
		color = builder.GenerateColor(project.Count.Experiments)
	}

	experiment, isNew, err := apiClient.CreateOrResumeExperiment(ctx, api.CreateExperimentParams{
		Username:    workspace,
		Project:     projectName,
		Name:        name,
		Description: rec.GetDescription(),
		Color:       color,
		JobType:     rec.GetJobType(),
		Group:       rec.GetGroup(),
		Tags:        rec.GetTags(),
		RunID:       rec.GetId(),
		Resume:      resumeMode,
	})
	if err != nil {
		return nil, fmt.Errorf("start: create or resume experiment: %w", err)
	}
	if experiment.CUID == "" || experiment.Name == "" {
		return nil, fmt.Errorf("%w: incomplete experiment metadata (cuid=%q name=%q)",
			ErrInvalidStartRecord, experiment.CUID, experiment.Name)
	}

	run.SetOnlineParams(runctx.OnlineParams{
		Username:       workspace,
		Project:        projectName,
		ProjectID:      project.CUID,
		ProjectVersion: int32(project.Version),
		ExperimentID:   experiment.CUID,
	})

	resume, err := restoreMetrics(ctx, apiClient, project.CUID, experiment, resumeMode, isNew)
	if err != nil {
		return nil, err
	}

	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.run = run
	s.metrics = resume.Metrics
	s.counter = counter.New(0)
	s.epoch = counter.New(resume.ConsoleEpoch)

	if err = s.open(sessionInit{
		cfg:          cfg,
		run:          run,
		http:         httpClient,
		api:          apiClient,
		workspace:    workspace,
		project:      projectName,
		projectID:    project.CUID,
		experimentID: experiment.CUID,
	}); err != nil {
		s.rollback()
		return nil, fmt.Errorf("start: assemble session: %w", err)
	}

	start := proto.Clone(rec).(*runv1.StartRecord)
	start.Name = name
	start.Color = color
	start.Workspace = workspace
	start.Project = projectName
	if err = s.commitStart(ctx, start); err != nil {
		s.rollback()
		return nil, err
	}

	sessionID, err := m.install(s)
	if err != nil {
		s.rollback()
		return nil, err
	}

	return &StartResult{
		SessionID:        sessionID,
		Run:              start,
		Path:             project.Path + "/" + firstNonEmpty(experiment.Slug, experiment.CUID),
		Name:             experiment.Name,
		GlobalStep:       resume.GlobalStep,
		GlobalSystemStep: resume.GlobalSystemStep,
		NewExperiment:    isNew,
	}, nil
}

// restoreMetrics 依据 resume 语义恢复指标状态；新实验或 never 模式返回空状态。
func restoreMetrics(
	ctx context.Context,
	client *api.Client,
	projectID string,
	experiment *api.Experiment,
	mode api.ResumeMode,
	isNew bool,
) (*metrics.ResumeResult, error) {
	empty := &metrics.ResumeResult{Metrics: metrics.New(), GlobalStep: -1, GlobalSystemStep: -1}
	if isNew || mode == api.ResumeNever {
		return empty, nil
	}
	summary, err := client.ExperimentSummary(ctx, projectID, experiment.CUID, experiment.CreatedAt)
	if err != nil {
		return nil, fmt.Errorf("start: fetch experiment summary: %w", err)
	}
	resume, err := metrics.NewFromSummary(summaryFromAPI(summary))
	if err != nil {
		return nil, fmt.Errorf("start: restore metrics: %w", err)
	}
	return resume, nil
}

// open 装配 store/intake/uploader/watcher/heartbeat；失败时由调用方执行 rollback。
func (s *Session) open(init sessionInit) error {
	if !init.cfg.SkipStore {
		runFile, err := init.run.RunFile()
		if err != nil {
			return fmt.Errorf("derive run file: %w", err)
		}
		if err = os.MkdirAll(init.run.RunDir, 0o755); err != nil {
			return fmt.Errorf("create run dir: %w", err)
		}
		writer, err := store.Create(runFile)
		if err != nil {
			return fmt.Errorf("create run file: %w", err)
		}
		s.runFile = runFile
		s.runFileCreated = true

		async, err := store.NewAsyncWriter(writer, store.AsyncOptions{
			Queue: queueSize(init.cfg.RecordBatch),
			OnError: func(err error) {
				s.storeDegraded.Store(true)
				console.Warningf("store degraded, records will be uploaded without a local copy: %v", err)
			},
		})
		if err != nil {
			_ = writer.Close()
			return fmt.Errorf("start store writer: %w", err)
		}
		s.store = async
	}

	s.tracker = uploader.NewTracker()
	b := budget.New(budget.Limits{})

	var storeIn chan<- store.Request
	var flushed func() uint64
	if s.store != nil {
		storeIn = s.store.In()
		flushed = func() uint64 {
			seq, _ := s.store.Flushed()
			return seq
		}
	}
	intake, err := buffer.New(buffer.Options{
		SkipStore:   init.cfg.SkipStore,
		Store:       storeIn,
		UploadQueue: queueSize(init.cfg.RecordBatch),
		Budget:      b,
		Degraded:    &s.storeDegraded,
		Flushed:     flushed,
		OnAdmit:     s.tracker.AddTotalRecords,
	})
	if err != nil {
		return fmt.Errorf("create intake: %w", err)
	}
	s.intake = intake

	senderOpts := sender.Options{
		API:          init.api,
		Username:     init.workspace,
		Project:      init.project,
		ProjectID:    init.projectID,
		ExperimentID: init.experimentID,
		Files: sender.FileOptions{
			MediaDir:  init.run.MediaDir(),
			FilesDir:  init.run.FilesDir(),
			SkipStore: init.cfg.SkipStore,
			SaveSize:  init.cfg.SaveSize,
			SaveSplit: init.cfg.SaveSplit,
			SavePart:  init.cfg.SavePart,
			SaveBatch: int(init.cfg.SaveBatch),
			Tracker:   s.tracker,
			Putter:    init.http,
		},
	}
	mediaSender, err := sender.NewMedia(senderOpts)
	if err != nil {
		return fmt.Errorf("create media sender: %w", err)
	}
	saveSender, err := sender.NewSave(senderOpts)
	if err != nil {
		return fmt.Errorf("create save sender: %w", err)
	}

	up, err := uploader.New(uploader.Options{
		In:     intake.UploadOut(),
		Budget: b,
		Senders: map[record.Type]uploader.Sender{
			record.TypeColumn: sender.NewColumn(senderOpts),
			record.TypeScalar: sender.NewScalar(senderOpts),
			record.TypeLog:    sender.NewLog(senderOpts),
			record.TypeMedia:  mediaSender,
			record.TypeSave:   saveSender,
		},
		Batch:    int(init.cfg.RecordBatch),
		Interval: time.Duration(init.cfg.RecordInterval * float64(time.Second)),
		Tracker:  s.tracker,
		Context:  s.ctx,
		OnDrained: func() {
			s.markDrained()
		},
	})
	if err != nil {
		return fmt.Errorf("create uploader: %w", err)
	}
	s.uploader = up
	s.drainFn = up.Wait
	s.releaseFn = func(ctx context.Context) error {
		return init.api.StopExperiment(ctx, api.StopExperimentParams{
			Username:     init.workspace,
			Project:      init.project,
			ExperimentID: init.experimentID,
			State:        s.stopState(),
		})
	}

	if init.cfg.SkipStore {
		s.watcher = watcher.NullWatcher{}
	} else {
		fileWatcher, err := watcher.New(init.run.FilesDir(), s.onFileChanged)
		if err != nil {
			return fmt.Errorf("create watcher: %w", err)
		}
		s.watcher = fileWatcher
	}

	s.tracker.SetState(operationv1.CoreState_CORE_STATE_RUNNING)
	s.heartbeat = newHeartbeat(init.experimentID, init.api.SendHeartbeat)
	s.heartbeat.start(s.ctx)
	if s.store != nil {
		// skip_store 不会有本地写失败，无需补传协程。
		s.startRecovery()
	}
	return nil
}

// commitStart 把 StartRecord 只写入本地 store；skip_store 时成为 no-op。
func (s *Session) commitStart(ctx context.Context, rec *runv1.StartRecord) error {
	ticket, err := s.intake.Reserve(ctx, buffer.Demand{Store: 1})
	if err != nil {
		return fmt.Errorf("start: reserve start record: %w", err)
	}
	if err := s.intake.Commit(ticket, []*recordv1.Record{builder.BuildStartRecord(rec)}); err != nil {
		return fmt.Errorf("start: commit start record: %w", err)
	}
	return nil
}

// onFileChanged 由 watcher 在 LIVE 文件稳定后调用，重新入队一条 save 记录。
// 防抖窗口内的残留由 watcher.Stop 在收尾时丢弃（状态门返回 ErrNotRunning）。
func (s *Session) onFileChanged(rec *savev1.SaveRecord) {
	if rec == nil {
		return
	}
	ctx := s.ctx
	if ctx == nil {
		ctx = context.Background()
	}
	err := s.commitPrepared(ctx, func() (preparedBatch, error) {
		built, buildErr := builder.BuildSaveRecord(s.counter, rec)
		if buildErr != nil {
			return preparedBatch{}, fmt.Errorf("%w: %v", ErrInvalidRecord, buildErr)
		}
		return preparedBatch{records: []*recordv1.Record{built}}, nil
	})
	if err == nil {
		return
	}
	if errors.Is(err, ErrNotRunning) {
		console.Debugf("finish: dropped pending file change for %q", rec.GetName())
		return
	}
	console.Warningf("watcher: re-enqueue %q: %v", rec.GetName(), err)
}

// rollback 释放装配期已创建的资源；所有字段都允许为 nil。
func (s *Session) rollback() {
	if s.cancel != nil {
		s.cancel()
	}
	if s.heartbeat != nil {
		s.heartbeat.stop()
	}
	if s.watcher != nil {
		s.watcher.Stop()
	}
	s.stopRecovery()
	if s.uploader != nil {
		_ = s.uploader.Close()
	}
	if s.intake != nil {
		ctx, cancel := context.WithTimeout(context.Background(), rollbackTimeout)
		_ = s.intake.Close(ctx)
		cancel()
	}
	if s.store != nil {
		_ = s.store.Close()
	}
	if s.runFileCreated {
		_ = os.Remove(s.runFile)
	}
}

// queueSize 是上传通道与 store 通道容量：单批记录最多因自动列翻倍。
func queueSize(batch int32) int { return queueFactor * int(batch) }

func firstNonEmpty(values ...string) string {
	for _, value := range values {
		if value != "" {
			return value
		}
	}
	return ""
}

// beginDrain 把会话从 RUNNING 切到 DRAINING；重复调用幂等。
// 已 FAILED 时返回原始失败，调用方不得把它当作成功。
func (s *Session) beginDrain() error {
	_, err := s.beginDrainOwned()
	return err
}

// beginDrainOwned 切换排空状态，并报告本次调用是否真正完成了 RUNNING → DRAINING。
// 只有赢得切换的调用方才能执行收尾序列，重复的 DeliverRunFinish 直接按幂等成功返回。
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

// finish 执行 DeliverRunFinish 的收尾序列：
//
//	DRAINING → watcher.Stop → 落盘 FinishRecord → 阻塞刷盘 → 补传 error_log →
//	Intake.Close（冲刷 END 暂存与 recovery）→ RequestFinish
//
// 本地刷盘不设超时，必须等 AsyncWriter 把已接收的记录全部写入磁盘；
// 刷盘失败只置 storeDegraded，由 intake 在 Close 时把未确认记录转入上传，
// 不让本地副本缺失阻断云端上传。FinishRecord 落盘或 Intake.Close 失败才置 FAILED。
func (s *Session) finish(ctx context.Context, rec *runv1.FinishRecord) error {
	started, err := s.beginDrainOwned()
	if err != nil {
		return err
	}
	if !started {
		return nil
	}

	// watcher.Stop 会等待在途回调并同步检查一次文件状态；此时状态已是 DRAINING，
	// 防抖窗口内的残留回调会被状态门拒绝并丢弃，与 Python 的 _store_finish 一致。
	if s.watcher != nil {
		s.watcher.Stop()
	}

	logRecord, commitErr := s.storeFinish(ctx, rec)
	if commitErr != nil {
		s.fail(commitErr)
	}

	if s.store != nil {
		if err := s.store.Close(); err != nil {
			s.storeDegraded.Store(true)
			console.Warningf("finish: local store did not flush completely: %v", err)
		}
	}

	// 回收队列中是更早的记录，error_log 不能插到它们前面。
	if logRecord != nil {
		if s.storeDegraded.Load() {
			s.drainRecovery(ctx)
		}
		if err := s.intake.CommitUploadOnly(ctx, logRecord); err != nil {
			console.Warningf("finish: upload error log: %v", err)
		}
	}

	s.stopRecovery()
	if err := s.closeIntake(ctx); err != nil {
		s.fail(err)
		return s.failure()
	}
	if commitErr != nil {
		return s.failure()
	}
	return nil
}

// closeIntake 冲刷 END 暂存与 recovery、关闭上传通道，并通知 uploader 收尾；
// 只有 Close 成功才允许 RequestFinish，避免 uploader 在通道未封闭时退出。
func (s *Session) closeIntake(ctx context.Context) error {
	if s.intake == nil {
		return nil
	}
	if err := s.intake.Close(ctx); err != nil {
		return fmt.Errorf("%w: close intake: %w", ErrRunFailed, err)
	}
	if s.uploader != nil {
		s.uploader.RequestFinish()
	}
	return nil
}

// storeFinish 构建并提交 FinishRecord（仅本地）；非 FINISHED 状态额外产出一条
// 崩溃 error_log，由调用方在刷盘完成后通过 CommitUploadOnly 补传。
func (s *Session) storeFinish(ctx context.Context, rec *runv1.FinishRecord) (*recordv1.Record, error) {
	finish := proto.Clone(rec).(*runv1.FinishRecord)

	s.mu.Lock()
	defer s.mu.Unlock()
	if RunState(s.state.Load()) != RunStateDraining {
		return nil, fmt.Errorf("%w: session left draining state", ErrNotRunning)
	}
	s.finishRecord = finish

	var logRecord *recordv1.Record
	if finish.GetState() != runv1.RunState_RUN_STATE_FINISHED {
		line := finish.GetError()
		if line == "" {
			line = finishUnknownError
		}
		logRecord = builder.BuildLogRecord(s.counter, s.epoch, &terminalv1.LogRecord{
			Line:      line,
			Level:     terminalv1.LogLevel_LOG_LEVEL_ERROR,
			Timestamp: finish.GetFinishedAt(),
		})
	}

	ticket, err := s.intake.Reserve(ctx, buffer.Demand{Store: 1})
	if err != nil {
		return nil, fmt.Errorf("reserve finish record: %w", err)
	}
	if err = s.intake.Commit(ticket, []*recordv1.Record{builder.BuildFinishRecord(finish)}); err != nil {
		return nil, fmt.Errorf("commit finish record: %w", err)
	}
	return logRecord, nil
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

// retryClose 在 Confirm 命中 FAILED 时重跑被中断的 Intake.Close，并等待 uploader
// 真正排空；成功后才把状态恢复为 DRAINED。Close 仍失败时保留 FAILED 供下次重试。
func (s *Session) retryClose(ctx context.Context) error {
	s.mu.Lock()
	failed := RunState(s.state.Load()) == RunStateFailed
	s.mu.Unlock()
	if !failed {
		return nil
	}
	if s.intake == nil {
		return s.failure()
	}
	if err := s.closeIntake(ctx); err != nil {
		s.fail(err)
		return s.failure()
	}
	if s.uploader != nil {
		if err := s.uploader.Wait(ctx); err != nil {
			return fmt.Errorf("confirm: wait for drain: %w", err)
		}
	}
	s.markRecovered()
	return nil
}

// stopState 返回上报给云端的实验终态；未记录 FinishRecord 时按正常结束处理。
func (s *Session) stopState() api.RunState {
	s.mu.Lock()
	defer s.mu.Unlock()
	switch s.finishRecord.GetState() {
	case runv1.RunState_RUN_STATE_CRASHED:
		return api.RunStateCrashed
	case runv1.RunState_RUN_STATE_ABORTED:
		return api.RunStateAborted
	default:
		return api.RunStateFinished
	}
}

// confirmMessage 是按本地副本完整性生成的 stop_experiment 失败文案。
func (s *Session) confirmMessage() string {
	switch {
	case s.cfg.SkipStore:
		return "Failed to finish run, and no local copy was kept (skip_store); the run state on the cloud is unconfirmed."
	case s.storeDegraded.Load():
		return "Failed to finish run, and the local copy is incomplete; the run state on the cloud is unconfirmed."
	default:
		return "Failed to finish run, but it has been saved locally."
	}
}

// waitDrained 在锁外等待会话排空，返回 nil 表示可以进入释放阶段。
func (s *Session) waitDrained(ctx context.Context) error {
	switch RunState(s.state.Load()) {
	case RunStateRunning:
		return ErrRunNotFinished
	case RunStateDrained:
		return nil
	case RunStateFailed:
		return s.failure()
	}

	if s.drainFn == nil {
		s.markDrained()
		return nil
	}
	if err := s.drainFn(ctx); err != nil {
		select {
		case <-s.failed:
			return s.failure()
		default:
			return err
		}
	}
	select {
	case <-s.drained:
		return nil
	case <-s.failed:
		return s.failure()
	case <-ctx.Done():
		return ctx.Err()
	}
}

// stateSnapshot 返回当前会话状态的并发安全快照。
func (s *Session) stateSnapshot() RunState {
	return RunState(s.state.Load())
}

// snapshot 返回上传进度快照；装配完成前 tracker 为空，返回空统计。
func (s *Session) snapshot() *operationv1.OperationStats {
	if s.tracker == nil {
		return &operationv1.OperationStats{}
	}
	return s.tracker.Snapshot()
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
