package syncer

import (
	"context"
	"errors"
	"fmt"
	"io"
	"math"
	"sync"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/runctx"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	"github.com/swanhubx/swanlab/core/internal/uploader/sender"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

const (
	// negCounterStart 与 Python pkg.builder.DEC_NUM 对齐：sync 合成的自动列与崩溃日志使用 <-100 编号，
	// 避免与 store 中的正数序号及 -3..-6 固定槽位冲突导致 recordBuffer 原地覆盖。
	negCounterStart        = -100
	defaultSyncBudgetBytes = 64 << 20
	queueFactor            = 2
	minQueueSize           = 16
	maxQueueSize           = 256
	maxBatchFloor          = 16

	finishUnknownError     = "run failed with unknown error while finish_record.error is not set"
	confirmOKMessage       = "OK"
	confirmFinishFailedMsg = "Failed to finish run, but all records have been uploaded to the server."
	confirmDrainFailedMsg  = "Failed to drain sync upload queue"
	missingStartRecordMsg  = "Failed to sync this run because its log file is missing the required start record. " +
		"The file may be corrupted, incomplete, or not a valid SwanLab run file."
	missingEmptyRunFileMsg = "Cannot scan records from %s"
	maxIntervalDuration    = time.Hour
)

type syncState int32

const (
	syncStateStarted syncState = iota
	syncStateFlushing
	syncStateRunning
	syncStateConfirming
	syncStateClosed
)

type syncSession struct {
	mu    sync.Mutex
	opMu  sync.Mutex
	state syncState

	cfg         config.Settings
	runFile     string
	startRecord *runv1.StartRecord
	workspace   string
	project     string
	runID       string
	filter      RangeFilter

	ctx           context.Context
	cancel        context.CancelFunc
	api           *api.Client
	run           *runctx.RunContext
	metrics       *metrics.RunMetrics
	epoch         *counter.Counter
	negCounter    *counter.Counter
	budget        *budget.Budget
	tracker       *uploader.Tracker
	uploader      *uploader.Uploader
	inCh          chan *recordv1.Record
	scanDone      chan struct{}
	newExperiment bool
	finishRecord  *runv1.FinishRecord
}

func newSyncSession(opts StartOptions, filter RangeFilter) (*syncSession, error) {
	if opts.Settings.RunDir == "" {
		return nil, fmt.Errorf("%w: %w", ErrInvalidSettings, config.ErrRunDirRequired)
	}
	runFile, err := runctx.FindSyncRunFile(opts.Settings.RunDir)
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
	}
	startRec, err := readFirstStartRecord(runFile)
	if err != nil {
		return nil, err
	}
	cfg, err := normalizeSyncSettings(opts.Settings, opts.RunID, startRec.GetId())
	if err != nil {
		return nil, err
	}

	return &syncSession{
		state:       syncStateStarted,
		cfg:         cfg,
		runFile:     runFile,
		startRecord: startRec,
		workspace:   opts.Workspace,
		project:     opts.Project,
		runID:       opts.RunID,
		filter:      filter,
	}, nil
}

func normalizeSyncSettings(cfg config.Settings, runIDOverride, fallbackRunID string) (config.Settings, error) {
	// sync 总是从本地读取并上传到云端，不受运行环境 offline/skip_store 配置影响。
	cfg.Mode = config.ModeOnline
	cfg.SkipStore = false
	cfg.RunID = firstNonEmpty(cfg.RunID, runIDOverride, fallbackRunID, "sync")
	if cfg.SaveSize <= 0 || cfg.SaveSplit <= 0 || cfg.SaveBatch <= 0 {
		return config.Settings{}, fmt.Errorf(
			"%w: save_size (%d), save_split (%d) and save_batch (%d) must be positive",
			ErrInvalidSettings, cfg.SaveSize, cfg.SaveSplit, cfg.SaveBatch,
		)
	}
	if err := cfg.Validate(); err != nil {
		return config.Settings{}, fmt.Errorf("%w: %w", ErrInvalidSettings, err)
	}
	return cfg, nil
}

func readFirstStartRecord(runFile string) (*runv1.StartRecord, error) {
	reader, err := store.Open(runFile)
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
	}
	defer func() { _ = reader.Close() }()

	payload, err := reader.Scan()
	if err != nil {
		if errors.Is(err, io.EOF) {
			return nil, fmt.Errorf("%w: "+missingEmptyRunFileMsg, ErrInvalidStartRecord, runFile)
		}
		return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
	}
	rec := &recordv1.Record{}
	if err := proto.Unmarshal(payload, rec); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalidStartRecord, err)
	}
	if rec.GetStart() == nil {
		return nil, fmt.Errorf("%w: %s", ErrInvalidStartRecord, missingStartRecordMsg)
	}
	return rec.GetStart(), nil
}

func (s *syncSession) flush(ctx context.Context, factory *client.Factory) (FlushResult, error) {
	s.opMu.Lock()
	defer s.opMu.Unlock()

	s.mu.Lock()
	if s.state != syncStateStarted {
		st := s.state
		s.mu.Unlock()
		return FlushResult{}, fmt.Errorf("%w: cannot flush from state %d", ErrInvalidState, st)
	}
	s.state = syncStateFlushing
	s.mu.Unlock()

	rollbackToStarted := func() {
		s.mu.Lock()
		if s.state == syncStateFlushing {
			s.state = syncStateStarted
		}
		s.mu.Unlock()
	}

	if s.cfg.APIHost == "" {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("%w: %w", ErrInvalidSettings, client.ErrMissingAPIHost)
	}
	if s.cfg.APIKey == "" {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("%w: %w", ErrInvalidSettings, client.ErrMissingAPIKey)
	}
	normalizedHost, err := client.NormalizeHost(s.cfg.APIHost)
	if err != nil {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: authenticate: %w", err)
	}

	cfg := s.cfg
	cfg.APIHost = normalizedHost
	run, err := runctx.New(cfg)
	if err != nil {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("%w: %w", ErrInvalidSettings, err)
	}

	httpClient, err := factory.Client(ctx, cfg.APIHost, cfg.APIKey)
	if err != nil {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: authenticate: %w", err)
	}
	apiClient := api.New(httpClient)

	start := proto.Clone(s.startRecord).(*runv1.StartRecord)
	if s.project != "" {
		start.Project = s.project
	}
	if s.workspace != "" {
		start.Workspace = s.workspace
	}
	if s.runID != "" {
		start.Id = s.runID
	}
	start.Resume = runv1.ResumeMode_RESUME_MODE_ALLOW
	start.Color = ""

	project, err := apiClient.GetOrCreateProject(ctx, start.GetWorkspace(), start.GetProject(), start.GetPublic())
	if err != nil {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: get project %q: %w", start.GetProject(), err)
	}
	workspace := firstNonEmpty(project.Group.Username, start.GetWorkspace(), apiClient.Username())
	projectName := firstNonEmpty(project.Name, start.GetProject())
	if workspace == "" || projectName == "" || project.CUID == "" || project.Path == "" {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf(
			"%w: incomplete project metadata (workspace=%q project=%q cuid=%q path=%q)",
			client.ErrInvalidResponse, workspace, projectName, project.CUID, project.Path,
		)
	}

	name := start.GetName()
	if name == "" {
		name = builder.GenerateName(project.Count.Experiments)
	}
	color := builder.GenerateColor(project.Count.Experiments)
	experiment, isNew, err := apiClient.CreateOrResumeExperiment(ctx, api.CreateExperimentParams{
		Username:    workspace,
		Project:     projectName,
		Name:        name,
		Description: start.GetDescription(),
		Color:       color,
		JobType:     start.GetJobType(),
		Group:       start.GetGroup(),
		Tags:        start.GetTags(),
		RunID:       start.GetId(),
		Resume:      api.ResumeAllow,
	})
	if err != nil {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: create or resume experiment: %w", err)
	}
	if experiment.CUID == "" || experiment.Name == "" {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("%w: incomplete experiment metadata (cuid=%q name=%q)",
			client.ErrInvalidResponse, experiment.CUID, experiment.Name)
	}

	run.SetOnlineParams(runctx.OnlineParams{
		Username:       workspace,
		Project:        projectName,
		ProjectID:      project.CUID,
		ProjectVersion: int32(project.Version),
		ExperimentID:   experiment.CUID,
	})

	var summary *metrics.Summary
	if !isNew {
		apiSummary, sumErr := apiClient.ExperimentSummary(ctx, project.CUID, experiment.CUID, experiment.CreatedAt)
		if sumErr != nil {
			rollbackToStarted()
			return FlushResult{}, fmt.Errorf("sync flush: fetch experiment summary: %w", sumErr)
		}
		summary = summaryFromAPI(apiSummary)
	}
	resumed, err := metrics.NewFromSummary(summary)
	if err != nil {
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: initialize metrics: %w", err)
	}
	if s.filter.Enabled() {
		resumed.Metrics.AllowReplay()
	}

	sessionCtx, cancel := context.WithCancel(context.Background())
	tracker := uploader.NewTracker()
	tracker.SetState(operationv1.CoreState_CORE_STATE_RUNNING)

	senderOpts := sender.Options{
		API:          apiClient,
		Username:     workspace,
		Project:      projectName,
		ProjectID:    project.CUID,
		ExperimentID: experiment.CUID,
		Files: sender.FileOptions{
			MediaDir:  run.MediaDir(),
			FilesDir:  run.FilesDir(),
			SkipStore: false,
			SaveSize:  cfg.SaveSize,
			SaveSplit: cfg.SaveSplit,
			SavePart:  cfg.SavePart,
			SaveBatch: int(cfg.SaveBatch),
			Tracker:   tracker,
			Putter:    httpClient,
		},
	}
	mediaSender, err := sender.NewMedia(senderOpts)
	if err != nil {
		cancel()
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: create media sender: %w", err)
	}
	saveSender, err := sender.NewSave(senderOpts)
	if err != nil {
		cancel()
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: create save sender: %w", err)
	}

	b := budget.New(budget.Limits{Bytes: defaultSyncBudgetBytes})
	inCh := make(chan *recordv1.Record, syncQueueSize(cfg.RecordBatch))
	up, err := uploader.New(uploader.Options{
		In:     inCh,
		Budget: b,
		Senders: map[record.Type]uploader.Sender{
			record.TypeColumn:    sender.NewColumn(senderOpts),
			record.TypeScalar:    sender.NewScalar(senderOpts),
			record.TypeLog:       sender.NewLog(senderOpts),
			record.TypeMedia:     mediaSender,
			record.TypeHistogram: sender.NewHistogram(senderOpts),
			record.TypeSave:      saveSender,
		},
		Batch:    int(cfg.RecordBatch),
		Interval: syncInterval(cfg.RecordInterval),
		Tracker:  tracker,
		Context:  sessionCtx,
		OnDrained: func() {
			tracker.SetState(operationv1.CoreState_CORE_STATE_FINISHED)
		},
	})
	if err != nil {
		cancel()
		rollbackToStarted()
		return FlushResult{}, fmt.Errorf("sync flush: create uploader: %w", err)
	}

	s.mu.Lock()
	if s.state != syncStateFlushing {
		s.mu.Unlock()
		cancel()
		_ = up.Close()
		return FlushResult{}, fmt.Errorf("%w: session closed during flush", ErrInvalidState)
	}
	s.ctx = sessionCtx
	s.cancel = cancel
	s.api = apiClient
	s.run = run
	s.metrics = resumed.Metrics
	s.epoch = counter.New(resumed.ConsoleEpoch)
	s.negCounter = counter.New(negCounterStart)
	s.budget = b
	s.tracker = tracker
	s.uploader = up
	s.inCh = inCh
	s.scanDone = make(chan struct{})
	s.newExperiment = isNew
	s.state = syncStateRunning
	s.mu.Unlock()

	go s.scanLoop()

	path := project.Path + "/" + firstNonEmpty(experiment.Slug, experiment.CUID)
	return FlushResult{Path: path}, nil
}

func (s *syncSession) stats() *operationv1.OperationStats {
	s.mu.Lock()
	tracker := s.tracker
	s.mu.Unlock()
	if tracker == nil {
		return &operationv1.OperationStats{State: operationv1.CoreState_CORE_STATE_NOT_STARTED}
	}
	return tracker.Snapshot()
}

func (s *syncSession) confirm(ctx context.Context) (ConfirmResult, error) {
	s.opMu.Lock()
	defer s.opMu.Unlock()

	s.mu.Lock()
	switch s.state {
	case syncStateRunning:
		s.state = syncStateConfirming
	case syncStateClosed:
		s.mu.Unlock()
		return ConfirmResult{Success: true, Message: confirmOKMessage}, nil
	default:
		st := s.state
		s.mu.Unlock()
		return ConfirmResult{}, fmt.Errorf("%w: cannot confirm from state %d", ErrInvalidState, st)
	}
	scanDone := s.scanDone
	up := s.uploader
	tracker := s.tracker
	s.mu.Unlock()

	select {
	case <-scanDone:
	case <-ctx.Done():
		s.restoreRunningOnCancel()
		return ConfirmResult{}, ctx.Err()
	}

	if err := up.Wait(ctx); err != nil {
		if ctx.Err() != nil {
			s.restoreRunningOnCancel()
			return ConfirmResult{}, ctx.Err()
		}
		console.Warningf("sync confirm: uploader wait failed: %v", err)
		s.cleanup()
		if tracker != nil {
			tracker.SetState(operationv1.CoreState_CORE_STATE_FINISHED)
		}
		return ConfirmResult{Success: false, Message: confirmDrainFailedMsg}, nil
	}

	var res ConfirmResult
	if state, ok := s.resolveFinishState(); !ok {
		res = ConfirmResult{Success: true, Message: confirmOKMessage}
	} else if err := s.stopExperiment(ctx, state); err != nil {
		console.Warningf("sync confirm: stop experiment: %v", err)
		res = ConfirmResult{Success: false, Message: confirmFinishFailedMsg}
	} else {
		res = ConfirmResult{Success: true, Message: confirmOKMessage}
	}

	s.cleanup()
	if tracker != nil {
		tracker.SetState(operationv1.CoreState_CORE_STATE_FINISHED)
	}
	return res, nil
}

func (s *syncSession) restoreRunningOnCancel() {
	s.mu.Lock()
	if s.state == syncStateConfirming {
		s.state = syncStateRunning
	}
	s.mu.Unlock()
}

func (s *syncSession) disconnect(ctx context.Context) error {
	s.mu.Lock()
	if s.state == syncStateClosed {
		s.mu.Unlock()
		return nil
	}
	s.state = syncStateClosed
	cancel := s.cancel
	scanDone := s.scanDone
	up := s.uploader
	tracker := s.tracker
	s.mu.Unlock()

	if cancel != nil {
		cancel()
	}
	var errs []error
	if scanDone != nil {
		select {
		case <-scanDone:
		case <-ctx.Done():
			errs = append(errs, ctx.Err())
		}
	}
	if up != nil {
		if err := up.Close(); err != nil && !errors.Is(err, uploader.ErrNotDrained) {
			errs = append(errs, err)
		}
	}
	if tracker != nil {
		tracker.SetState(operationv1.CoreState_CORE_STATE_FINISHED)
	}
	return errors.Join(errs...)
}

func (s *syncSession) cleanup() {
	s.mu.Lock()
	s.state = syncStateClosed
	cancel := s.cancel
	up := s.uploader
	s.mu.Unlock()

	if cancel != nil {
		cancel()
	}
	if up != nil {
		_ = up.Close()
	}
}

func (s *syncSession) resolveFinishState() (api.RunState, bool) {
	if s.finishRecord != nil {
		return runStateFromProto(s.finishRecord.GetState()), true
	}
	if s.newExperiment {
		return api.RunStateCrashed, true
	}
	return "", false
}

func (s *syncSession) stopExperiment(ctx context.Context, state api.RunState) error {
	online, ok := s.run.Online()
	if !ok {
		return errors.New("sync session missing online params")
	}
	return s.api.StopExperiment(ctx, api.StopExperimentParams{
		Username:     online.Username,
		Project:      online.Project,
		ExperimentID: online.ExperimentID,
		State:        state,
	})
}

func runStateFromProto(state runv1.RunState) api.RunState {
	switch state {
	case runv1.RunState_RUN_STATE_CRASHED:
		return api.RunStateCrashed
	case runv1.RunState_RUN_STATE_ABORTED:
		return api.RunStateAborted
	default:
		return api.RunStateFinished
	}
}

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

func syncQueueSize(recordBatch int32) int {
	size := max(int(recordBatch), maxBatchFloor) * queueFactor
	size = max(size, minQueueSize)
	return min(size, maxQueueSize)
}

func syncInterval(seconds float64) time.Duration {
	if math.IsInf(seconds, 0) || seconds <= 0 {
		return maxIntervalDuration
	}
	ns := seconds * float64(time.Second)
	if ns > float64(maxIntervalDuration) {
		return maxIntervalDuration
	}
	return max(time.Duration(ns), time.Millisecond)
}

func firstNonEmpty(values ...string) string {
	for _, v := range values {
		if v != "" {
			return v
		}
	}
	return ""
}
