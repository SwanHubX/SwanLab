// 本文件承载 DeliverRunStart 的云端编排、指标恢复与启动记录提交。
package manager

import (
	"context"
	"errors"
	"fmt"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/buffer"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

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

// summaryFromAPI 把云端摘要转换为指标恢复层使用的纯 Go 结构。
// 两侧字段一一对应，nil 输入保持 nil，交给 metrics.NewFromSummary 处理。
func summaryFromAPI(s *api.Summary) *metrics.Summary {
	if s == nil {
		return nil
	}
	return &metrics.Summary{
		Log:    summaryItemsFromAPI(s.Log),
		Media:  summaryItemsFromAPI(s.Media),
		Scalar: summaryItemsFromAPI(s.Scalar),
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
