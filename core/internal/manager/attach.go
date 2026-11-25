// 本文件负责在后端不可用时先记录到本地，并在后端恢复后创建云端实验、上传历史记录。
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
	// 连续认证失败达到这个次数后，认为凭据或权限无效，停止重试。
	attachAuthFailLimit = 3

	// 连接重试的默认等待时间和单次请求超时。
	attachInitialBackoff = 500 * time.Millisecond
	attachMaxBackoff     = 30 * time.Second
	attachJitterMax      = 500 * time.Millisecond
	attachAttemptTimeout = 30 * time.Second
	// 站点状态通知意外关闭时，改为按这个间隔检查主服务是否恢复。
	attachServerPoll = time.Second
)

// errAttachStopped 表示会话已经结束或正在退出，不能再保存云端资源。
// 这是正常的停止条件，不应记录为云端连接永久失败。
var errAttachStopped = errors.New("cloud attach stopped by session lifecycle")

// attachTiming 控制连接重试的等待时间和超时；零值字段使用默认值。
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

// attachParams 保存后台创建云端实验和上传资源所需的参数。
type attachParams struct {
	factory    *client.Factory
	start      *runv1.StartRecord
	resumeMode api.ResumeMode
	timing     attachTiming
	// 名称和颜色在第一次拿到项目信息后确定，后续重试复用，避免被重新生成。
	name  string
	color string
}

// startAttach 为会话启动后台任务，负责创建云端实验并保存上传资源。
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

// stopAttach 取消并等待后台连接任务退出；未启动或已停止时可以重复调用。
func (s *Session) stopAttach() {
	s.cloudMu.Lock()
	cancel := s.attachCancel
	s.cloudMu.Unlock()
	if cancel != nil {
		cancel()
	}
	s.attachWG.Wait()
}

// waitAttach 等待后台云端连接结束。没有启动连接或已经连接成功时立即返回。
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
			// 连接任务因会话结束被取消：这是主动停止，不是云端永久失败。
			return fmt.Errorf("%w: cloud attach stopped before completion", context.Canceled)
		}
	case <-ctx.Done():
		return ctx.Err()
	}
}

// runAttach 按退避重试创建云端实验，直到成功、遇到不可恢复错误或任务被取消。
// 站点健康检查恢复只用于提前唤醒；多个会话同时恢复时随机错开请求。
func (s *Session) runAttach(ctx context.Context, params attachParams) {
	// 退出时关闭 attachDone，让等待连接结果的结束流程可以继续；
	// 不修改 cloudState：既不能标记为成功，也不能误报为永久失败。
	defer s.closeAttachDone()
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

		retryable := false
		if client.IsAuthenticationError(err) {
			authFailures++
			if authFailures >= attachAuthFailLimit {
				console.Errorf("cloud attach failed permanently after %d authentication errors: %v", authFailures, err)
				s.markCloudFailed(err)
				return
			}
			retryable = true
		} else if isRetryableAttachError(err) {
			authFailures = 0
			retryable = true
		}

		if !retryable {
			// 非认证且不可重试的错误（如参数无效、resume=must 找不到实验、响应格式异常）不再重试。
			console.Errorf("cloud attach failed permanently: %v", err)
			s.markCloudFailed(err)
			return
		}

		s.noteAttachError(err)
		delay := nextAttachDelay(err, backoff, timing.maxBackoff)
		if !waitAttachBackoff(ctx, delay, siteUpdates, s.isServerDegraded, timing.jitterMax) {
			return
		}
		backoff = min(backoff*2, timing.maxBackoff)
	}
}

// tryAttach 执行一次完整连接：认证、获取或创建项目、创建实验、创建上传资源并保存到会话。
// 只有会话仍在运行时才启动后台上传；DRAINING/FAILED 由结束流程统一处理积压记录。
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
	if err = requireProjectMetadata(workspace, projectName, project.CUID, project.Path); err != nil {
		return err
	}

	if ctx.Err() != nil {
		return ctx.Err()
	}
	// 会话已经完成结束时不再创建实验，避免留下没有会话管理的云端实验。
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

	// 站点读取侧异常且本地存储可用时，先暂停上传，避免向异常的服务发请求。
	if s.shouldSuspendAttachedUploader() {
		res.uploader.Suspend()
	}

	// 在 Session.mu 内同时检查状态、保存云端资源并决定是否启动后台上传：
	// - RUNNING：保存资源并启动后台上传；
	// - DRAINING/FAILED：只保存资源，积压记录由结束流程统一上传；
	// - 其他状态或 ctx 已取消：关闭刚创建的资源，返回 errAttachStopped。
	// FAILED 也必须保存资源，否则已经创建的实验无法在 Confirm 中停止。
	// 连接任务不获取 finishMu，避免和等待连接结果的结束流程互相等待。
	var startedBackfill bool
	stopped := false
	s.mu.Lock()
	state := RunState(s.state.Load())
	switch {
	case ctx.Err() != nil:
		stopped = true
	case state == RunStateRunning, state == RunStateDraining, state == RunStateFailed:
		s.run.SetOnlineParams(runctx.OnlineParams{
			Username:       workspace,
			Project:        projectName,
			ProjectID:      project.CUID,
			ProjectVersion: int32(project.Version),
			ExperimentID:   experiment.CUID,
		})
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
	return s.canSuspendUploads() && s.isSiteReadDegraded()
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

// waitForServerRecovery 等待主服务恢复。恢复后随机等待一小段时间，避免多个会话同时请求。
func waitForServerRecovery(
	ctx context.Context,
	siteUpdates <-chan struct{},
	serverDegraded func() bool,
	jitterMax time.Duration,
) bool {
	for serverDegraded() {
		if siteUpdates == nil {
			// 没有站点通知时按固定间隔检查，保证连接任务不会就此退出。
			timer := time.NewTimer(attachServerPoll)
			select {
			case <-ctx.Done():
				timer.Stop()
				return false
			case <-timer.C:
			}
			continue
		}
		select {
		case <-ctx.Done():
			return false
		case _, ok := <-siteUpdates:
			if !ok {
				// 订阅意外关闭后改为定时检查，继续等待主服务恢复。
				siteUpdates = nil
			}
		}
	}
	return waitJitter(ctx, jitterMax)
}

// waitAttachBackoff 等待下一次重试；收到服务恢复通知时提前唤醒并随机错开。
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
				// 订阅通道关闭后只按定时器等待，避免空转。
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
