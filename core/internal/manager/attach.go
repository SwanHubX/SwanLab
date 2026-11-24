// 本文件实现初始化阶段后端不可达时的异步云端接入与历史补传启动逻辑。
package manager

import (
	"context"
	"errors"
	"fmt"
	"math/rand/v2"
	"time"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/cron/health"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

const (
	// 连续认证失败次数上限，达到后认定凭据或权限失效，停止重试。
	attachAuthFailLimit = 3

	// 接入退避与错峰抖动的默认时间参数。
	attachInitialBackoff = 500 * time.Millisecond
	attachMaxBackoff     = 30 * time.Second
	attachJitterMax      = 500 * time.Millisecond
	attachAttemptTimeout = 30 * time.Second
	// 站点订阅通道意外关闭后，等待主服务恢复的轮询间隔。
	attachServerPoll = time.Second
)

// errAttachStopped 表示会话生命周期已不允许发布云端资源（例如已进入 DRAINED 或正在退出），
// 属于正常停止而不是可上报的永久接入失败。
var errAttachStopped = errors.New("cloud attach stopped by session lifecycle")

// attachTiming 控制接入重试的退避与超时参数；零值字段使用生产默认常量。
type attachTiming struct {
	initialBackoff time.Duration
	maxBackoff     time.Duration
	jitterMax      time.Duration
	attemptTimeout time.Duration
}

func (t attachTiming) withDefaults() attachTiming {
	if t.initialBackoff <= 0 {
		t.initialBackoff = attachInitialBackoff
	}
	if t.maxBackoff <= 0 {
		t.maxBackoff = attachMaxBackoff
	}
	if t.jitterMax <= 0 {
		t.jitterMax = attachJitterMax
	}
	if t.attemptTimeout <= 0 {
		t.attemptTimeout = attachAttemptTimeout
	}
	return t
}

// attachParams 保存后台创建云端实验并挂载上传资源所需的全部静态参数。
type attachParams struct {
	factory    *client.Factory
	start      *runv1.StartRecord
	resumeMode api.ResumeMode
	timing     attachTiming
	// name/color 在首次成功读取项目统计后固定，重试时不再随项目实验数变化重新生成。
	name  string
	color string
}

// startAttach 启动单个会话的后台接入协程；受会话生命周期与 stopAttach 控制。
func (s *Session) startAttach(params attachParams) {
	if params.name == "" {
		params.name = params.start.GetName()
	}
	if params.color == "" {
		params.color = params.start.GetColor()
	}
	parent := s.ctx
	if parent == nil {
		parent = context.Background()
	}
	attachCtx, cancel := context.WithCancel(parent)
	s.cloudMu.Lock()
	s.attachCancel = cancel
	s.cloudMu.Unlock()

	s.attachWG.Add(1)
	go func() {
		defer s.attachWG.Done()
		s.runAttach(attachCtx, params)
	}()
}

// stopAttach 取消并等待后台接入协程退出；未启动或已停止时可重复调用。
func (s *Session) stopAttach() {
	s.cloudMu.Lock()
	cancel := s.attachCancel
	s.cloudMu.Unlock()
	if cancel != nil {
		cancel()
	}
	s.attachWG.Wait()
}

// waitAttach 等待延迟接入完成或永久失败；未启用延迟接入或已接入时立即返回。
func (s *Session) waitAttach(ctx context.Context) error {
	s.cloudMu.Lock()
	done := s.attachDone
	s.cloudMu.Unlock()
	if done == nil {
		return nil
	}
	select {
	case <-done:
		state, _, err := s.cloudStatus()
		switch state {
		case operationv1.CloudState_CLOUD_STATE_FAILED:
			return err
		case operationv1.CloudState_CLOUD_STATE_ATTACHED:
			return nil
		default:
			// 接入协程因退出/收尾取消而停止：这是上下文级中断，不是云端永久失败。
			return fmt.Errorf("%w: cloud attach stopped before completion", context.Canceled)
		}
	case <-ctx.Done():
		return ctx.Err()
	}
}

// runAttach 周期性尝试接入云端，直到成功、遇到不可恢复错误或会话取消。
// 重试由自身退避驱动；站点健康探针恢复时仅作为立即唤醒信号，多会话恢复时加入随机抖动错峰。
func (s *Session) runAttach(ctx context.Context, params attachParams) {
	// 无论以何种方式退出，都要释放 waitAttach，避免收尾协程永久等待接入结果。
	defer s.markAttachStopped()
	timing := params.timing.withDefaults()
	var siteUpdates <-chan struct{}
	if s.site != nil {
		var unsubscribe func()
		siteUpdates, unsubscribe = s.site.Subscribe()
		defer unsubscribe()
	}

	backoff := timing.initialBackoff
	authFailures := 0

	for {
		if ctx.Err() != nil {
			return
		}
		if s.isServerDegraded() {
			if !waitForServerRecovery(ctx, siteUpdates, s.isServerDegraded, timing.jitterMax) {
				return
			}
			backoff = timing.initialBackoff
		}

		err := s.tryAttach(ctx, &params, timing.attemptTimeout)
		if err == nil {
			return
		}
		if errors.Is(err, errAttachStopped) || ctx.Err() != nil {
			return
		}

		if client.IsAuthenticationError(err) {
			authFailures++
			if authFailures >= attachAuthFailLimit {
				console.Errorf("cloud attach failed permanently after %d authentication errors: %v", authFailures, err)
				s.markCloudFailed(err)
				return
			}
			s.noteAttachError(err)
			delay := nextAttachDelay(err, backoff, timing.maxBackoff)
			if !waitAttachBackoff(ctx, delay, siteUpdates, s.isServerDegraded, timing.jitterMax) {
				return
			}
			backoff = min(backoff*2, timing.maxBackoff)
			continue
		}

		if isRetryableAttachError(err) {
			authFailures = 0
			s.noteAttachError(err)
			delay := nextAttachDelay(err, backoff, timing.maxBackoff)
			if !waitAttachBackoff(ctx, delay, siteUpdates, s.isServerDegraded, timing.jitterMax) {
				return
			}
			backoff = min(backoff*2, timing.maxBackoff)
			continue
		}

		// 非认证且不可重试的错误（如参数拒绝、resume=must 找不到实验、响应格式异常）不再重试。
		console.Errorf("cloud attach failed permanently: %v", err)
		s.markCloudFailed(err)
		return
	}
}

// tryAttach 执行一次有序接入：认证 -> 建项目/实验 -> 装配云端资源 -> 原子挂载 -> 条件启动补传。
func (s *Session) tryAttach(ctx context.Context, params *attachParams, attemptTimeout time.Duration) error {
	attemptCtx, cancel := context.WithTimeout(ctx, attemptTimeout)
	defer cancel()

	httpClient, err := params.factory.Client(attemptCtx, s.cfg.APIHost, s.cfg.APIKey)
	if err != nil {
		return fmt.Errorf("attach: authenticate: %w", err)
	}
	apiClient := api.New(httpClient)

	project, err := apiClient.GetOrCreateProject(
		attemptCtx,
		params.start.GetWorkspace(),
		params.start.GetProject(),
		params.start.GetPublic(),
	)
	if err != nil {
		return fmt.Errorf("attach: get project %q: %w", params.start.GetProject(), err)
	}
	workspace := firstNonEmpty(project.Group.Username, params.start.GetWorkspace(), apiClient.Username())
	projectName := firstNonEmpty(project.Name, params.start.GetProject())
	if workspace == "" || projectName == "" || project.CUID == "" || project.Path == "" {
		return fmt.Errorf(
			"%w: incomplete project metadata (workspace=%q project=%q cuid=%q path=%q)",
			client.ErrInvalidResponse, workspace, projectName, project.CUID, project.Path,
		)
	}

	if ctx.Err() != nil {
		return ctx.Err()
	}
	// 已 DRAINED 的会话不再创建实验，避免退出路径留下无人接管的云端身份。
	if s.stateSnapshot() == RunStateDrained {
		return errAttachStopped
	}
	if params.name == "" {
		params.name = builder.GenerateName(project.Count.Experiments)
	}
	if params.color == "" {
		params.color = builder.GenerateColor(project.Count.Experiments)
	}
	runID := firstNonEmpty(params.start.GetId(), s.cfg.RunID)

	experiment, _, err := apiClient.CreateOrResumeExperiment(attemptCtx, api.CreateExperimentParams{
		Username:    workspace,
		Project:     projectName,
		Name:        params.name,
		Description: params.start.GetDescription(),
		Color:       params.color,
		JobType:     params.start.GetJobType(),
		Group:       params.start.GetGroup(),
		Tags:        params.start.GetTags(),
		RunID:       runID,
		Resume:      params.resumeMode,
	})
	if err != nil {
		return fmt.Errorf("attach: create or resume experiment: %w", err)
	}
	if experiment.CUID == "" || experiment.Name == "" {
		return fmt.Errorf(
			"%w: incomplete experiment metadata (cuid=%q name=%q)",
			client.ErrInvalidResponse, experiment.CUID, experiment.Name,
		)
	}

	if ctx.Err() != nil {
		return ctx.Err()
	}

	s.run.SetOnlineParams(runctx.OnlineParams{
		Username:       workspace,
		Project:        projectName,
		ProjectID:      project.CUID,
		ProjectVersion: int32(project.Version),
		ExperimentID:   experiment.CUID,
	})

	res, err := s.buildCloudResources(cloudAttachParams{
		http:         httpClient,
		api:          apiClient,
		workspace:    workspace,
		project:      projectName,
		projectID:    project.CUID,
		experimentID: experiment.CUID,
	})
	if err != nil {
		return fmt.Errorf("attach: build resources: %w", err)
	}

	// 站点非写入侧降级且本地存储可用时，挂载前保持暂停，避免向异常的存储服务发请求。
	if s.shouldSuspendAttachedUploader() {
		res.uploader.Suspend()
	}

	// 状态检查、挂载云端资源与“仅在运行中启动补传”共用 Session.mu，避免与收尾停补传竞态；
	// 接入协程绝不获取 finishMu，防止与等待接入的收尾协程死锁。
	//
	// FAILED 也允许挂载：DeliverRunFinish 可能在实验创建期间失败，此时 Confirm 会通过
	// resumeFinish 统一补传和排空，不能让已创建的实验失去 releaseFn。只有退出路径
	// （ctx 已取消或已进入 DRAINED）才放弃挂载，并返回可识别的停止错误。
	var startedBackfill bool
	stopped := false
	s.mu.Lock()
	state := RunState(s.state.Load())
	switch {
	case ctx.Err() != nil:
		stopped = true
	case state == RunStateRunning, state == RunStateDraining, state == RunStateFailed:
		s.setCloudResources(res)
		if state == RunStateRunning {
			s.backfill.Start()
			startedBackfill = true
		}
	default:
		stopped = true
	}
	s.mu.Unlock()
	if stopped {
		_ = res.uploader.Close()
		return errAttachStopped
	}

	cloudPath := project.Path + "/" + firstNonEmpty(experiment.Slug, experiment.CUID)
	res.heartbeat.Start(s.ctx)
	s.markCloudAttached(cloudPath)
	if startedBackfill {
		s.backfill.Wake()
	}
	s.wakeHealth()
	console.Infof("cloud attach succeeded: experiment=%s path=%s", experiment.CUID, cloudPath)
	return nil
}

func (s *Session) isServerDegraded() bool {
	if s.site == nil {
		return false
	}
	snap := s.site.Snapshot()
	return snap.State == health.StateDegraded && snap.Reason == health.ReasonProbeServer
}

func (s *Session) shouldSuspendAttachedUploader() bool {
	if s.site == nil || s.cfg.SkipStore || s.storeDegraded.Load() {
		return false
	}
	snap := s.site.Snapshot()
	return snap.State == health.StateDegraded && snap.Reason != health.ReasonWrite
}

func isRetryableAttachError(err error) bool {
	return client.IsRetryable(err) || errors.Is(err, context.DeadlineExceeded)
}

func nextAttachDelay(err error, backoff, maxBackoff time.Duration) time.Duration {
	delay := min(backoff, maxBackoff)
	var apiErr *client.APIError
	if errors.As(err, &apiErr) && apiErr.RetryAfter > delay {
		delay = min(apiErr.RetryAfter, maxBackoff)
	}
	return delay
}

// waitForServerRecovery 等待站点主服务探针脱离降级，并在唤醒后加入随机抖动错峰。
func waitForServerRecovery(
	ctx context.Context,
	siteUpdates <-chan struct{},
	serverDegraded func() bool,
	jitterMax time.Duration,
) bool {
	if siteUpdates == nil {
		// 无订阅时退化为定时轮询，保证接入协程不会静默退出。
		for serverDegraded() {
			timer := time.NewTimer(attachServerPoll)
			select {
			case <-ctx.Done():
				timer.Stop()
				return false
			case <-timer.C:
			}
		}
		return waitJitter(ctx, jitterMax)
	}
	for serverDegraded() {
		select {
		case <-ctx.Done():
			return false
		case _, ok := <-siteUpdates:
			if !ok {
				// 订阅通道意外关闭：退化为定时轮询，继续等待主服务恢复。
				return waitForServerRecovery(ctx, nil, serverDegraded, jitterMax)
			}
		}
	}
	return waitJitter(ctx, jitterMax)
}

// waitAttachBackoff 按退避时长等待下一次重试；若期间收到服务端恢复通知则立即唤醒并错峰。
func waitAttachBackoff(
	ctx context.Context,
	delay time.Duration,
	siteUpdates <-chan struct{},
	serverDegraded func() bool,
	jitterMax time.Duration,
) bool {
	if delay <= 0 {
		return ctx.Err() == nil
	}
	timer := time.NewTimer(delay)
	defer timer.Stop()
	for {
		select {
		case <-ctx.Done():
			return false
		case <-timer.C:
			return true
		case _, ok := <-siteUpdates:
			if !ok {
				// 订阅通道意外关闭时降级为纯定时器等待，避免忙轮询。
				siteUpdates = nil
				continue
			}
			if !serverDegraded() {
				return waitJitter(ctx, jitterMax)
			}
		}
	}
}

func waitJitter(ctx context.Context, maxJitter time.Duration) bool {
	if maxJitter <= 0 {
		return ctx.Err() == nil
	}
	jitter := rand.N(maxJitter)
	if jitter <= 0 {
		return ctx.Err() == nil
	}
	timer := time.NewTimer(jitter)
	defer timer.Stop()
	select {
	case <-ctx.Done():
		return false
	case <-timer.C:
		return true
	}
}
