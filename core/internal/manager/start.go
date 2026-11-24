// 本文件承载 DeliverRunStart 的云端编排、指标恢复与启动记录提交。
package manager

import (
	"context"
	"errors"
	"fmt"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/cron/health"
	"github.com/swanhubx/swanlab/core/internal/manager/probe"
	"github.com/swanhubx/swanlab/core/internal/manager/resume"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

type cloudBootstrap struct {
	apiHost    string
	http       *client.Client
	api        *api.Client
	workspace  string
	project    string
	projectID  string
	projectVer int32
	experiment *api.Experiment
	path       string
	name       string
	color      string
	isNew      bool
	healthOpts health.SiteOptions
	resumeRes  *metrics.ResumeResult
}

func buildHealthOpts(
	cfg config.Settings,
	serverProbe, houseProbe health.ProbeFunc,
	serverDegraded, houseDegraded bool,
) health.SiteOptions {
	// 0 值节奏与阈值由 health 补默认。
	return health.SiteOptions{
		Key:                   probe.SiteKey(cfg.APIHost, cfg.APIKey),
		Interval:              probe.SecondsDuration(cfg.HealthProbeInterval),
		Timeout:               probe.SecondsDuration(cfg.HealthProbeTimeout),
		DegradeThreshold:      int(cfg.HealthDegradeThreshold),
		RecoverThreshold:      int(cfg.HealthRecoverThreshold),
		Server:                serverProbe,
		House:                 houseProbe,
		InitialServerDegraded: serverDegraded,
		InitialHouseDegraded:  houseDegraded,
	}
}

func (m *Manager) bootstrapCloud(
	ctx context.Context,
	cfg config.Settings,
	rec *runv1.StartRecord,
	resumeMode api.ResumeMode,
	local *metrics.RebuildResult,
) (*cloudBootstrap, error) {
	normalizedHost, err := client.NormalizeHost(cfg.APIHost)
	if err != nil {
		return nil, fmt.Errorf("start: authenticate: %w", err)
	}
	cfg.APIHost = normalizedHost

	// 预检与认证客户端使用相同的规范化 host；server 不可达时不再探测 house。
	serverProbe := probe.Site(cfg.APIHost, probe.ServerPath)
	houseProbe := probe.Site(cfg.APIHost, probe.HousePath)
	serverOK, houseOK := probe.Startup(ctx, serverProbe, houseProbe)
	if !serverOK {
		return nil, fmt.Errorf("%w: %s", ErrServerUnreachable, probe.URL(cfg.APIHost, probe.ServerPath))
	}
	healthOpts := buildHealthOpts(cfg, serverProbe, houseProbe, false, !houseOK)

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
			client.ErrInvalidResponse, workspace, projectName, project.CUID, project.Path)
	}
	if err = resume.ValidateStart(local, resume.StartIdentity{
		RunID:               firstNonEmpty(rec.GetId(), cfg.RunID),
		Workspace:           workspace,
		Project:             projectName,
		AllowEmptyWorkspace: true,
	}); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
	}

	name := rec.GetName()
	color := rec.GetColor()
	if local != nil && local.Start != nil {
		// 续传时本地 StartRecord 是既有实验身份，优先于本次请求里的默认值。
		if local.Start.GetName() != "" {
			name = local.Start.GetName()
		}
		if local.Start.GetColor() != "" {
			color = local.Start.GetColor()
		}
	}
	if name == "" {
		name = builder.GenerateName(project.Count.Experiments)
	}
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
			client.ErrInvalidResponse, experiment.CUID, experiment.Name)
	}

	resumeRes, err := restoreMetrics(ctx, apiClient, project.CUID, experiment, resumeMode, isNew, local)
	if err != nil {
		return nil, err
	}

	return &cloudBootstrap{
		apiHost:    normalizedHost,
		http:       httpClient,
		api:        apiClient,
		workspace:  workspace,
		project:    projectName,
		projectID:  project.CUID,
		projectVer: int32(project.Version),
		experiment: experiment,
		path:       project.Path + "/" + firstNonEmpty(experiment.Slug, experiment.CUID),
		name:       name,
		color:      color,
		isNew:      isNew,
		healthOpts: healthOpts,
		resumeRes:  resumeRes,
	}, nil
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
	local, err := resume.LoadLocalRebuild(cfg, run)
	if err != nil {
		if errors.Is(err, resume.ErrInvalidRunFile) {
			return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
		}
		return nil, err
	}
	if cfg.NormalizedMode() != config.ModeOnline {
		return m.deliverLocalStart(ctx, cfg, rec, run, local)
	}
	resumeMode, err := resumeModeFromProto(rec.GetResume())
	if err != nil {
		return nil, fmt.Errorf("%w: %w", ErrInvalidStartRecord, err)
	}
	if cfg.RecoverOnline {
		return m.deliverRecoverStart(ctx, cfg, rec, run, resumeMode, local)
	}
	boot, err := m.bootstrapCloud(ctx, cfg, rec, resumeMode, local)
	if err != nil {
		return nil, err
	}
	cfg.APIHost = boot.apiHost

	run.SetOnlineParams(runctx.OnlineParams{
		Username:       boot.workspace,
		Project:        boot.project,
		ProjectID:      boot.projectID,
		ProjectVersion: boot.projectVer,
		ExperimentID:   boot.experiment.CUID,
	})

	var startCounter int64
	if local != nil {
		startCounter = local.Counter
	}
	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.run = run
	s.metrics = boot.resumeRes.Metrics
	s.counter = counter.New(startCounter)
	s.epoch = counter.New(boot.resumeRes.ConsoleEpoch)

	if err = s.open(sessionInit{
		cfg:          cfg,
		run:          run,
		http:         boot.http,
		api:          boot.api,
		workspace:    boot.workspace,
		project:      boot.project,
		projectID:    boot.projectID,
		experimentID: boot.experiment.CUID,
		healthReg:    m.healthReg,
		healthOpts:   boot.healthOpts,
		resume:       local,
	}); err != nil {
		s.rollback()
		return nil, fmt.Errorf("start: assemble session: %w", err)
	}

	start := proto.Clone(rec).(*runv1.StartRecord)
	start.Name = boot.name
	start.Color = boot.color
	start.Workspace = boot.workspace
	start.Project = boot.project
	if local == nil || local.Start == nil {
		if err = s.commitStart(ctx, start); err != nil {
			s.rollback()
			return nil, err
		}
	}

	sessionID, err := m.install(s)
	if err != nil {
		s.rollback()
		return nil, err
	}

	return &StartResult{
		SessionID:        sessionID,
		Run:              start,
		Path:             boot.path,
		Name:             boot.experiment.Name,
		GlobalStep:       boot.resumeRes.GlobalStep,
		GlobalSystemStep: boot.resumeRes.GlobalSystemStep,
		NewExperiment:    boot.isNew,
	}, nil
}

// deliverRecoverStart 在云端暂不可用时先装配本地在线会话，并启动后台异步接入任务。
func (m *Manager) deliverRecoverStart(
	ctx context.Context,
	cfg config.Settings,
	rec *runv1.StartRecord,
	run *runctx.RunContext,
	resumeMode api.ResumeMode,
	local *metrics.RebuildResult,
) (*StartResult, error) {
	if resumeMode != api.ResumeNever && (local == nil || local.Records == 0) {
		return nil, fmt.Errorf("%w: recover_online with resume=%s requires an existing local run file",
			ErrInvalidStartRecord, resumeMode)
	}
	if err := resume.ValidateStart(local, resume.StartIdentity{
		RunID:               firstNonEmpty(rec.GetId(), cfg.RunID),
		Workspace:           rec.GetWorkspace(),
		Project:             rec.GetProject(),
		AllowEmptyWorkspace: true,
	}); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
	}
	normalizedHost, err := client.NormalizeHost(cfg.APIHost)
	if err != nil {
		return nil, fmt.Errorf("start: authenticate: %w", err)
	}
	cfg.APIHost = normalizedHost

	serverProbe := probe.Site(cfg.APIHost, probe.ServerPath)
	houseProbe := probe.Site(cfg.APIHost, probe.HousePath)
	serverOK, houseOK := probe.Startup(ctx, serverProbe, houseProbe)
	healthOpts := buildHealthOpts(cfg, serverProbe, houseProbe, !serverOK, serverOK && !houseOK)

	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.run = run
	s.initDeferredAttach()

	globalStep, globalSystemStep := int64(-1), int64(-1)
	newExperiment := true
	if local != nil {
		s.metrics = local.Metrics
		s.counter = counter.New(local.Counter)
		s.epoch = counter.New(local.Epoch)
		globalStep = local.GlobalStep
		globalSystemStep = local.GlobalSystemStep
		newExperiment = false
	} else {
		s.metrics = metrics.New()
		s.counter = counter.New(0)
		s.epoch = counter.New(0)
	}

	if err = s.open(sessionInit{
		cfg:        cfg,
		run:        run,
		healthReg:  m.healthReg,
		healthOpts: healthOpts,
		resume:     local,
	}); err != nil {
		s.rollback()
		return nil, fmt.Errorf("start: assemble recover session: %w", err)
	}

	start := proto.Clone(rec).(*runv1.StartRecord)
	if local != nil && local.Start != nil {
		// 本地已有未完成的 StartRecord 时以它为准，保证恢复前后实验身份稳定。
		if local.Start.GetName() != "" {
			start.Name = local.Start.GetName()
		}
		if local.Start.GetColor() != "" {
			start.Color = local.Start.GetColor()
		}
	}
	if local == nil || local.Start == nil {
		if err = s.commitStart(ctx, start); err != nil {
			s.rollback()
			return nil, err
		}
	}

	// 注册表容量只能在 install 内最终确认；先启动接入会在安装失败时留下无人接管的
	// 云端实验。install 成功后再启动后台接入，rollback 就不会再遇到已创建的云端身份。
	sessionID, err := m.install(s)
	if err != nil {
		s.rollback()
		return nil, err
	}
	s.startAttach(attachParams{
		factory:    m.factory,
		start:      proto.Clone(start).(*runv1.StartRecord),
		resumeMode: resumeMode,
		timing:     m.attachTiming,
	})

	return &StartResult{
		SessionID:        sessionID,
		Run:              start,
		Name:             start.GetName(),
		GlobalStep:       globalStep,
		GlobalSystemStep: globalSystemStep,
		NewExperiment:    newExperiment,
	}, nil
}

// deliverLocalStart 为 local/offline 装配仅落盘的会话：不认证、不访问云端，
// 不创建 uploader/heartbeat，所有记录只进入本地 store。
func (m *Manager) deliverLocalStart(
	ctx context.Context,
	cfg config.Settings,
	rec *runv1.StartRecord,
	run *runctx.RunContext,
	local *metrics.RebuildResult,
) (*StartResult, error) {
	if err := resume.ValidateStart(local, resume.StartIdentity{
		RunID:     firstNonEmpty(rec.GetId(), cfg.RunID),
		Workspace: rec.GetWorkspace(),
		Project:   rec.GetProject(),
	}); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
	}
	s := newSession(sessionOptions{})
	s.ctx, s.cancel = context.WithCancel(context.Background())
	s.cfg = cfg
	s.run = run
	globalStep, globalSystemStep := int64(-1), int64(-1)
	newExperiment := true
	if local != nil {
		s.metrics = local.Metrics
		s.counter = counter.New(local.Counter)
		s.epoch = counter.New(local.Epoch)
		globalStep = local.GlobalStep
		globalSystemStep = local.GlobalSystemStep
		newExperiment = false
	} else {
		s.metrics = metrics.New()
		s.counter = counter.New(0)
		s.epoch = counter.New(0)
	}

	if err := s.open(sessionInit{cfg: cfg, run: run, resume: local}); err != nil {
		s.rollback()
		return nil, fmt.Errorf("start: assemble local session: %w", err)
	}

	start := proto.Clone(rec).(*runv1.StartRecord)
	if local == nil || local.Start == nil {
		if err := s.commitStart(ctx, start); err != nil {
			s.rollback()
			return nil, err
		}
	}
	sessionID, err := m.install(s)
	if err != nil {
		s.rollback()
		return nil, err
	}
	return &StartResult{
		SessionID:        sessionID,
		Run:              start,
		GlobalStep:       globalStep,
		GlobalSystemStep: globalSystemStep,
		NewExperiment:    newExperiment,
	}, nil
}

// restoreMetrics 依据 resume 语义合并本地重建结果与云端摘要；两者均无时返回空状态。
func restoreMetrics(
	ctx context.Context,
	client *api.Client,
	projectID string,
	experiment *api.Experiment,
	mode api.ResumeMode,
	isNew bool,
	local *metrics.RebuildResult,
) (*metrics.ResumeResult, error) {
	var cloud *metrics.Summary
	if !isNew && mode != api.ResumeNever {
		summary, err := client.ExperimentSummary(ctx, projectID, experiment.CUID, experiment.CreatedAt)
		if err != nil {
			return nil, fmt.Errorf("start: fetch experiment summary: %w", err)
		}
		cloud = summaryFromAPI(summary)
	}
	if local != nil {
		res, err := metrics.MergeSummary(local, cloud)
		if err != nil {
			return nil, fmt.Errorf("start: merge metrics: %w", err)
		}
		return res, nil
	}
	res, err := metrics.NewFromSummary(cloud)
	if err != nil {
		return nil, fmt.Errorf("start: restore metrics: %w", err)
	}
	return res, nil
}

// summaryFromAPI 把云端摘要转换为指标恢复层使用的纯 Go 结构。
// 两侧字段一一对应，nil 输入保持 nil，交给 metrics.NewFromSummary 处理。
func summaryFromAPI(s *api.Summary) *metrics.Summary {
	if s == nil {
		return nil
	}
	return &metrics.Summary{
		Log:       summaryItemsFromAPI(s.Log),
		Media:     summaryItemsFromAPI(s.Media),
		Scalar:    summaryItemsFromAPI(s.Scalar),
		Histogram: summaryItemsFromAPI(s.Histogram),
	}
}

func summaryItemsFromAPI(items []api.ColumnStep) []metrics.SummaryItem {
	if len(items) == 0 {
		return nil
	}
	out := make([]metrics.SummaryItem, len(items))
	for i, item := range items {
		out[i] = metrics.SummaryItem{Key: item.Key, Step: item.Step}
	}
	return out
}

// resumeModeFromProto 把 proto 枚举转换为云端接口使用的恢复策略。
// 未知取值属于协议错误，不能静默降级为 allow/never。
func resumeModeFromProto(mode runv1.ResumeMode) (api.ResumeMode, error) {
	switch mode {
	case runv1.ResumeMode_RESUME_MODE_NEVER:
		return api.ResumeNever, nil
	case runv1.ResumeMode_RESUME_MODE_ALLOW:
		return api.ResumeAllow, nil
	case runv1.ResumeMode_RESUME_MODE_MUST:
		return api.ResumeMust, nil
	default:
		return "", fmt.Errorf("unknown resume mode %d", mode)
	}
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
