package service

import (
	"context"
	"fmt"
	"net"
	"strings"
	"sync"
	"testing"
	"time"

	"go.uber.org/goleak"
	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/grpc/stats"
	"google.golang.org/grpc/status"
	"google.golang.org/grpc/test/bufconn"
	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/manager"
	"github.com/swanhubx/swanlab/core/internal/server"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	settingsv1 "github.com/swanhubx/swanlab/core/proto/swanlab/settings/core/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

const (
	bufconnSize      = 1 << 20
	callTimeout      = 2 * time.Second
	noShutdownWindow = 200 * time.Millisecond
)

func TestMain(m *testing.M) {
	goleak.VerifyTestMain(m)
}

// fakeManager 只实现 service 用到的编排接口，用于验证 RPC 契约与错误映射。
type fakeManager struct {
	mu       sync.Mutex
	next     int
	sessions map[string]manager.RunState
	starts   []fakeStart

	maxSessions    int
	startErr       error
	forceEmptyPath bool
	finishErr      error
	stateErr       error
	confirmErr     error
	upsertErr      error
	disconnectErr  error

	confirmResult manager.ConfirmResult

	upserts      []fakeUpsert
	finishes     []fakeFinish
	disconnected []string

	onStart          func()
	onDisconnect     func(context.Context, string)
	disconnectNotify chan string
}

type fakeFinish struct {
	sessionID string
	record    *runv1.FinishRecord
}

type fakeUpsert struct {
	kind      string
	sessionID string
	count     int
}

type fakeStart struct {
	settings config.Settings
	record   *runv1.StartRecord
}

func newFakeManager() *fakeManager {
	return &fakeManager{
		sessions:      make(map[string]manager.RunState),
		confirmResult: manager.ConfirmResult{Message: "OK"},
	}
}

func (f *fakeManager) DeliverRunStart(
	_ context.Context,
	settings config.Settings,
	record *runv1.StartRecord,
) (*manager.StartResult, error) {
	f.mu.Lock()
	hook := f.onStart
	f.mu.Unlock()
	if hook != nil {
		hook()
	}

	f.mu.Lock()
	defer f.mu.Unlock()
	f.starts = append(f.starts, fakeStart{settings: settings, record: record})
	if f.startErr != nil {
		return nil, f.startErr
	}
	if f.maxSessions > 0 && len(f.sessions) >= f.maxSessions {
		return nil, fmt.Errorf("%w: limit %d", manager.ErrTooManySessions, f.maxSessions)
	}
	f.next++
	sessionID := fmt.Sprintf("session-%d", f.next)
	f.sessions[sessionID] = manager.RunStateRunning
	path := "/alice/demo/run-1"
	if settings.NormalizedMode() != config.ModeOnline || settings.RecoverOnline || f.forceEmptyPath {
		path = ""
	}
	return &manager.StartResult{
		SessionID:     sessionID,
		Run:           record,
		Path:          path,
		Name:          "demo",
		GlobalStep:    -1,
		NewExperiment: true,
	}, nil
}

func (f *fakeManager) recordUpsert(kind, sessionID string, count int) error {
	f.mu.Lock()
	defer f.mu.Unlock()
	f.upserts = append(f.upserts, fakeUpsert{kind: kind, sessionID: sessionID, count: count})
	if f.upsertErr != nil {
		return f.upsertErr
	}
	if sessionID == "" {
		return manager.ErrEmptySessionID
	}
	if _, ok := f.sessions[sessionID]; !ok {
		return manager.ErrSessionNotFound
	}
	return nil
}

func (f *fakeManager) UpsertColumns(_ context.Context, sessionID string, records []*columnv1.ColumnRecord) error {
	return f.recordUpsert("columns", sessionID, len(records))
}

func (f *fakeManager) UpsertScalars(_ context.Context, sessionID string, records []*datav1.ScalarRecord) error {
	return f.recordUpsert("scalars", sessionID, len(records))
}

func (f *fakeManager) UpsertMedia(_ context.Context, sessionID string, records []*datav1.MediaRecord) error {
	return f.recordUpsert("media", sessionID, len(records))
}

func (f *fakeManager) UpsertHistograms(_ context.Context, sessionID string, records []*datav1.HistogramRecord) error {
	return f.recordUpsert("histograms", sessionID, len(records))
}

func (f *fakeManager) UpsertLogs(_ context.Context, sessionID string, records []*terminalv1.LogRecord) error {
	return f.recordUpsert("logs", sessionID, len(records))
}

func (f *fakeManager) UpsertSaves(_ context.Context, sessionID string, records []*savev1.SaveRecord) error {
	return f.recordUpsert("saves", sessionID, len(records))
}

func (f *fakeManager) DeliverRunFinish(_ context.Context, sessionID string, record *runv1.FinishRecord) error {
	f.mu.Lock()
	defer f.mu.Unlock()
	f.finishes = append(f.finishes, fakeFinish{sessionID: sessionID, record: record})
	if record == nil {
		return manager.ErrInvalidFinishRecord
	}
	if f.finishErr != nil {
		return f.finishErr
	}
	if sessionID == "" {
		return manager.ErrEmptySessionID
	}
	state, ok := f.sessions[sessionID]
	if !ok {
		return manager.ErrSessionNotFound
	}
	if state == manager.RunStateFailed {
		return manager.ErrRunFailed
	}
	f.sessions[sessionID] = manager.RunStateDraining
	return nil
}

func (f *fakeManager) State(sessionID string) (manager.RunState, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.stateErr != nil {
		return 0, f.stateErr
	}
	if sessionID == "" {
		return 0, manager.ErrEmptySessionID
	}
	state, ok := f.sessions[sessionID]
	if !ok {
		return 0, manager.ErrSessionNotFound
	}
	return state, nil
}

func (f *fakeManager) Stats(sessionID string) (*operationv1.OperationStats, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.stateErr != nil {
		return nil, f.stateErr
	}
	if sessionID == "" {
		return nil, manager.ErrEmptySessionID
	}
	if _, ok := f.sessions[sessionID]; !ok {
		return nil, manager.ErrSessionNotFound
	}
	return &operationv1.OperationStats{TotalRecords: 7, Rate: 1.5}, nil
}

func (f *fakeManager) Confirm(_ context.Context, sessionID string) (manager.ConfirmResult, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.confirmErr != nil {
		return manager.ConfirmResult{}, f.confirmErr
	}
	if sessionID == "" {
		return manager.ConfirmResult{}, manager.ErrEmptySessionID
	}
	state, ok := f.sessions[sessionID]
	if !ok {
		return manager.ConfirmResult{}, manager.ErrSessionNotFound
	}
	switch state {
	case manager.RunStateRunning:
		return manager.ConfirmResult{}, manager.ErrRunNotFinished
	case manager.RunStateFailed:
		return manager.ConfirmResult{}, manager.ErrRunFailed
	}
	if f.confirmResult.ReleaseFailed {
		return f.confirmResult, nil
	}
	delete(f.sessions, sessionID)
	return f.confirmResult, nil
}

func (f *fakeManager) DisconnectSession(ctx context.Context, sessionID string) error {
	f.mu.Lock()
	hook := f.onDisconnect
	notify := f.disconnectNotify
	f.mu.Unlock()
	if hook != nil {
		hook(ctx, sessionID)
	}

	f.mu.Lock()
	f.disconnected = append(f.disconnected, sessionID)
	delete(f.sessions, sessionID)
	err := f.disconnectErr
	f.mu.Unlock()

	if notify != nil {
		notify <- sessionID
	}
	return err
}

func (f *fakeManager) disconnectedSnapshot() []string {
	f.mu.Lock()
	defer f.mu.Unlock()
	return append([]string(nil), f.disconnected...)
}

type testEnv struct {
	client corev1.CoreServiceClient
	ctrl   *server.Controller
	runs   *fakeManager
}

func dialBufconn(t *testing.T, lis *bufconn.Listener) (*grpc.ClientConn, corev1.CoreServiceClient) {
	t.Helper()
	conn, err := grpc.NewClient(
		"passthrough:///bufnet",
		grpc.WithTransportCredentials(insecure.NewCredentials()),
		grpc.WithContextDialer(func(context.Context, string) (net.Conn, error) {
			return lis.DialContext(context.Background())
		}),
	)
	if err != nil {
		t.Fatalf("grpc.NewClient: %v", err)
	}
	t.Cleanup(func() { _ = conn.Close() })
	return conn, corev1.NewCoreServiceClient(conn)
}

func newTestEnv(t *testing.T) *testEnv {
	t.Helper()
	g := grpc.NewServer()
	ctrl := server.NewController(g, 2*time.Second)
	runs := newFakeManager()
	NewCoreService(ctrl, runs).Register(g)
	lis := bufconn.Listen(bufconnSize)
	go func() { _ = g.Serve(lis) }()
	t.Cleanup(func() {
		ctrl.Shutdown("test-cleanup")
		<-ctrl.Done()
	})
	_, client := dialBufconn(t, lis)
	return &testEnv{client: client, ctrl: ctrl, runs: runs}
}

func callCtx(t *testing.T) context.Context {
	t.Helper()
	ctx, cancel := context.WithTimeout(context.Background(), callTimeout)
	t.Cleanup(cancel)
	return ctx
}

func spinup(t *testing.T, env *testEnv) {
	t.Helper()
	if _, err := env.client.SpinupService(callCtx(t), &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
}

func validSettings(t *testing.T) *settingsv1.CoreSettings {
	t.Helper()
	return &settingsv1.CoreSettings{
		RunId:          "run-1",
		RunDir:         t.TempDir(),
		RecordBatch:    1000,
		RecordInterval: 1.5,
		SaveSize:       50 << 30,
		SaveSplit:      100 << 20,
		SavePart:       32 << 20,
		SaveBatch:      100,
		ApiKey:         "smoke-key",
		ApiHost:        "http://127.0.0.1:1",
	}
}

func startRun(t *testing.T, env *testEnv, record *runv1.StartRecord) string {
	t.Helper()
	if record == nil {
		record = &runv1.StartRecord{}
	}
	resp, err := env.client.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  record,
	})
	if err != nil || !resp.GetSuccess() || resp.GetSessionId() == "" {
		t.Fatalf("DeliverRunStart: resp=%v err=%v", resp, err)
	}
	return resp.GetSessionId()
}

func TestSpinupAndTeardown(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	spinup(t, env)
	if _, err := env.client.TeardownService(callCtx(t), &corev1.TeardownServiceRequest{}); err != nil {
		t.Fatalf("TeardownService: %v", err)
	}
	select {
	case <-env.ctrl.Done():
	case <-time.After(callTimeout):
		t.Fatal("shutdown not completed after teardown")
	}
}

func TestRunRPCsRejectedBeforeReady(t *testing.T) {
	env := newTestEnv(t)
	ctx := callCtx(t)
	for name, call := range map[string]func() error{
		"start": func() error {
			_, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
				CoreSettings: validSettings(t),
				StartRecord:  &runv1.StartRecord{},
			})
			return err
		},
		"upsert": func() error {
			_, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{SessionId: "any"})
			return err
		},
		"stats": func() error {
			_, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: "any"})
			return err
		},
		"finish": func() error {
			_, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{SessionId: "any"})
			return err
		},
		"confirm": func() error {
			_, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{SessionId: "any"})
			return err
		},
	} {
		if err := call(); status.Code(err) != codes.FailedPrecondition {
			t.Fatalf("%s before ready: err=%v, want FailedPrecondition", name, err)
		}
	}
}

func TestDeliverRunStartLocalAllowsEmptyPath(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)

	settings := validSettings(t)
	settings.Mode = settingsv1.CoreMode_CORE_MODE_LOCAL
	settings.ApiKey = ""
	settings.ApiHost = ""
	resp, err := env.client.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{
		CoreSettings: settings,
		StartRecord:  &runv1.StartRecord{Project: "demo", Id: "run-1"},
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	if !resp.GetSuccess() || resp.GetSessionId() == "" {
		t.Fatalf("local start response = %v", resp)
	}
	if len(env.runs.starts) != 1 || env.runs.starts[0].settings.Mode != config.ModeLocal {
		t.Fatalf("manager settings = %+v", env.runs.starts)
	}
}

func TestDeliverRunStartRecoverOnlineAllowsEmptyPath(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)

	settings := validSettings(t)
	settings.SkipStore = false
	settings.RecoverOnline = true

	resp, err := env.client.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{
		CoreSettings: settings,
		StartRecord:  &runv1.StartRecord{Project: "demo", Id: "run-recover"},
	})
	if err != nil {
		t.Fatalf("DeliverRunStart with RecoverOnline: %v", err)
	}
	if !resp.GetSuccess() || resp.GetSessionId() == "" || resp.GetPath() != "" {
		t.Fatalf("recover start response = %v, want success with empty path", resp)
	}

	// 未开启 RecoverOnline 时，online 模式返回空 Path 仍应被拒绝。
	env.runs.forceEmptyPath = true
	noRecover := validSettings(t)
	noRecover.SkipStore = false
	if _, err = env.client.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{
		CoreSettings: noRecover,
		StartRecord:  &runv1.StartRecord{Project: "demo", Id: "run-normal"},
	}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("normal online empty path err=%v, want InvalidArgument", err)
	}
}

func TestDeliverRunStartValidatesRequest(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)

	if _, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
	}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("missing start_record err=%v, want InvalidArgument", err)
	}
	if _, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		StartRecord: &runv1.StartRecord{},
	}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("missing core_settings err=%v, want InvalidArgument", err)
	}

	settings := validSettings(t)
	settings.SaveSize = 0
	if _, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: settings,
		StartRecord:  &runv1.StartRecord{},
	}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("zero save_size err=%v, want InvalidArgument", err)
	}
}

func TestDeliverRunStartPassesSettings(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	startRun(t, env, &runv1.StartRecord{Project: "demo"})
	if len(env.runs.starts) != 1 {
		t.Fatalf("start calls = %d, want 1", len(env.runs.starts))
	}
	got := env.runs.starts[0]
	if got.settings.RunID != "run-1" || got.settings.RecordBatch != 1000 {
		t.Fatalf("settings not forwarded: %+v", got.settings)
	}
	if got.record.GetProject() != "demo" {
		t.Fatalf("start record not forwarded: %v", got.record)
	}
}

func TestDeliverRunStartMapsCloudErrors(t *testing.T) {
	cases := []struct {
		name string
		err  error
		want codes.Code
	}{
		{"auth", client.ErrNotAuthenticated, codes.Unauthenticated},
		{"not found", &client.APIError{StatusCode: 404}, codes.NotFound},
		{"unavailable", &client.APIError{StatusCode: 503, Retryable: true}, codes.Unavailable},
		{"invalid start", fmt.Errorf("start: %w", manager.ErrInvalidStartRecord), codes.InvalidArgument},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			env := newTestEnv(t)
			spinup(t, env)
			env.runs.startErr = tc.err
			_, err := env.client.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{
				CoreSettings: validSettings(t),
				StartRecord:  &runv1.StartRecord{},
			})
			if code := status.Code(err); code != tc.want {
				t.Fatalf("code = %v, want %v (err=%v)", code, tc.want, err)
			}
		})
	}
}

func TestSettingsFromProto(t *testing.T) {
	in := &settingsv1.CoreSettings{
		RunId:          "run-1",
		RunDir:         "/tmp/run",
		SectionRule:    2,
		RecordBatch:    500,
		RecordInterval: 0.5,
		SaveSize:       1,
		SaveSplit:      2,
		SavePart:       3,
		SaveBatch:      4,
		SkipStore:      true,
		ApiKey:         "key",
		ApiHost:        "http://127.0.0.1:1",
		WebHost:        "http://127.0.0.1:2",

		UploadBudgetRecords:    0,
		UploadBudgetBytes:      0,
		HealthProbeInterval:    1.5,
		HealthProbeTimeout:     0.75,
		HealthDegradeThreshold: 5,
		HealthRecoverThreshold: 4,
	}
	got, err := settingsFromProto(in)
	if err != nil {
		t.Fatalf("settingsFromProto: %v", err)
	}
	want := config.Settings{
		RunID:          "run-1",
		RunDir:         "/tmp/run",
		Mode:           config.ModeOnline,
		SectionRule:    2,
		RecordBatch:    500,
		RecordInterval: 0.5,
		SaveSize:       1,
		SaveSplit:      2,
		SavePart:       3,
		SaveBatch:      4,
		SkipStore:      true,
		APIKey:         "key",
		APIHost:        "http://127.0.0.1:1",
		WebHost:        "http://127.0.0.1:2",

		UploadBudgetRecords:    0,
		UploadBudgetBytes:      0,
		HealthProbeInterval:    1.5,
		HealthProbeTimeout:     0.75,
		HealthDegradeThreshold: 5,
		HealthRecoverThreshold: 4,
	}
	if got != want {
		t.Fatalf("settings = %+v, want %+v", got, want)
	}

	local := &settingsv1.CoreSettings{
		RunId:          "run-1",
		RunDir:         "/tmp/run",
		Mode:           settingsv1.CoreMode_CORE_MODE_LOCAL,
		RecordBatch:    500,
		RecordInterval: 0.5,
		SaveSize:       1,
		SaveSplit:      2,
		SavePart:       3,
		SaveBatch:      4,
	}
	gotLocal, err := settingsFromProto(local)
	if err != nil {
		t.Fatalf("local settings must not require cloud credentials: %v", err)
	}
	if gotLocal.Mode != config.ModeLocal {
		t.Fatalf("mode = %q, want %q", gotLocal.Mode, config.ModeLocal)
	}

	if _, nilErr := settingsFromProto(nil); nilErr == nil {
		t.Fatal("nil settings must be rejected")
	}

	cases := []struct {
		name   string
		mutate func(*settingsv1.CoreSettings)
	}{
		{"missing api key", func(s *settingsv1.CoreSettings) { s.ApiKey = "" }},
		{"missing api host", func(s *settingsv1.CoreSettings) { s.ApiHost = "" }},
		{"zero save size", func(s *settingsv1.CoreSettings) { s.SaveSize = 0 }},
		{"zero save batch", func(s *settingsv1.CoreSettings) { s.SaveBatch = 0 }},
		{"zero record batch", func(s *settingsv1.CoreSettings) { s.RecordBatch = 0 }},
		{"zero save part", func(s *settingsv1.CoreSettings) { s.SavePart = 0 }},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			clone := proto.Clone(in).(*settingsv1.CoreSettings)
			tc.mutate(clone)
			if _, caseErr := settingsFromProto(clone); caseErr == nil {
				t.Fatal("expected an error")
			}
		})
	}

	budgeted := proto.Clone(in).(*settingsv1.CoreSettings)
	budgeted.UploadBudgetRecords = 50
	budgeted.UploadBudgetBytes = 1024
	budgetSettings, err := settingsFromProto(budgeted)
	if err != nil {
		t.Fatalf("non-zero upload budget must be accepted: %v", err)
	}
	if budgetSettings.UploadBudgetRecords != 50 || budgetSettings.UploadBudgetBytes != 1024 {
		t.Fatalf("budget = %d/%d, want 50/1024",
			budgetSettings.UploadBudgetRecords, budgetSettings.UploadBudgetBytes)
	}

	recoverOnline := proto.Clone(in).(*settingsv1.CoreSettings)
	recoverOnline.SkipStore = false
	recoverOnline.RecoverOnline = true
	gotRecover, err := settingsFromProto(recoverOnline)
	if err != nil {
		t.Fatalf("recover_online in online mode without skip_store must be accepted: %v", err)
	}
	if !gotRecover.RecoverOnline {
		t.Fatal("RecoverOnline = false, want true")
	}

	recoverRejectCases := []struct {
		name   string
		mutate func(*settingsv1.CoreSettings)
	}{
		{"recover_online with skip_store", func(s *settingsv1.CoreSettings) { s.SkipStore = true; s.RecoverOnline = true }},
		{"recover_online with local", func(s *settingsv1.CoreSettings) {
			s.Mode = settingsv1.CoreMode_CORE_MODE_LOCAL
			s.SkipStore = false
			s.RecoverOnline = true
		}},
		{"recover_online with offline", func(s *settingsv1.CoreSettings) {
			s.Mode = settingsv1.CoreMode_CORE_MODE_OFFLINE
			s.SkipStore = false
			s.RecoverOnline = true
		}},
		{"recover_online with disabled", func(s *settingsv1.CoreSettings) {
			s.Mode = settingsv1.CoreMode_CORE_MODE_DISABLED
			s.SkipStore = false
			s.RecoverOnline = true
		}},
	}
	for _, tc := range recoverRejectCases {
		t.Run(tc.name, func(t *testing.T) {
			clone := proto.Clone(in).(*settingsv1.CoreSettings)
			tc.mutate(clone)
			_, caseErr := settingsFromProto(clone)
			if caseErr == nil || !strings.Contains(caseErr.Error(), "recover_online") {
				t.Fatalf("settingsFromProto() = %v, want error mentioning recover_online", caseErr)
			}
		})
	}
}

func TestUpsertRPCsForwardRecords(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env, &runv1.StartRecord{})
	ctx := callCtx(t)

	calls := []struct {
		kind string
		call func() error
	}{
		{"columns", func() error {
			_, err := env.client.UpsertColumns(ctx, &corev1.UpsertColumnsRequest{
				SessionId: handle,
				Columns:   []*columnv1.ColumnRecord{{ColumnKey: "loss"}},
			})
			return err
		}},
		{"scalars", func() error {
			_, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{
				SessionId: handle,
				Data:      []*datav1.ScalarRecord{{Key: "loss"}},
			})
			return err
		}},
		{"media", func() error {
			_, err := env.client.UpsertMedia(ctx, &corev1.UpsertMediaRequest{
				SessionId: handle,
				Data:      []*datav1.MediaRecord{{Key: "image"}},
			})
			return err
		}},
		{"histograms", func() error {
			_, err := env.client.UpsertHistograms(ctx, &corev1.UpsertHistogramsRequest{
				SessionId: handle,
				Data:      []*datav1.HistogramRecord{{Key: "weights"}},
			})
			return err
		}},
		{"logs", func() error {
			_, err := env.client.UpsertLogs(ctx, &corev1.UpsertLogsRequest{
				SessionId: handle,
				Logs:      []*terminalv1.LogRecord{{Line: "hello"}},
			})
			return err
		}},
		{"saves", func() error {
			_, err := env.client.UpsertSaves(ctx, &corev1.UpsertSavesRequest{
				SessionId: handle,
				Saves:     []*savev1.SaveRecord{{Name: "model.pt"}},
			})
			return err
		}},
	}
	for _, tc := range calls {
		if err := tc.call(); err != nil {
			t.Fatalf("%s: %v", tc.kind, err)
		}
	}
	if len(env.runs.upserts) != len(calls) {
		t.Fatalf("forwarded upserts = %d, want %d", len(env.runs.upserts), len(calls))
	}
	for i, tc := range calls {
		got := env.runs.upserts[i]
		if got.kind != tc.kind || got.sessionID != handle || got.count != 1 {
			t.Fatalf("upsert %d = %+v, want kind=%s session=%s count=1", i, got, tc.kind, handle)
		}
	}
}

func TestUpsertStoreUnavailableMapsFailedPrecondition(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env, &runv1.StartRecord{})

	env.runs.mu.Lock()
	env.runs.upsertErr = manager.ErrStoreUnavailable
	env.runs.mu.Unlock()

	_, err := env.client.UpsertScalars(callCtx(t), &corev1.UpsertScalarsRequest{
		SessionId: handle,
		Data:      []*datav1.ScalarRecord{{Key: "loss"}},
	})
	if code := status.Code(err); code != codes.FailedPrecondition {
		t.Fatalf("code = %v, want FailedPrecondition (err=%v)", code, err)
	}
}

func TestUpsertEmptyListsAreNoops(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env, &runv1.StartRecord{})
	ctx := callCtx(t)

	// gRPC 不区分「缺省」与「空 repeated」，两者都应作为 no-op 成功。
	if _, err := env.client.UpsertColumns(ctx, &corev1.UpsertColumnsRequest{SessionId: handle}); err != nil {
		t.Fatalf("empty UpsertColumns: %v", err)
	}
	if _, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{SessionId: handle}); err != nil {
		t.Fatalf("empty UpsertScalars: %v", err)
	}
	if _, err := env.client.UpsertMedia(ctx, &corev1.UpsertMediaRequest{SessionId: handle}); err != nil {
		t.Fatalf("empty UpsertMedia: %v", err)
	}
	if _, err := env.client.UpsertHistograms(ctx, &corev1.UpsertHistogramsRequest{SessionId: handle}); err != nil {
		t.Fatalf("empty UpsertHistograms: %v", err)
	}
	if _, err := env.client.UpsertLogs(ctx, &corev1.UpsertLogsRequest{SessionId: handle}); err != nil {
		t.Fatalf("empty UpsertLogs: %v", err)
	}
	if _, err := env.client.UpsertSaves(ctx, &corev1.UpsertSavesRequest{SessionId: handle}); err != nil {
		t.Fatalf("empty UpsertSaves: %v", err)
	}
}

func TestUpsertRPCsMapErrors(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	handle := startRun(t, env, &runv1.StartRecord{})
	ctx := callCtx(t)

	tests := []struct {
		name string
		err  error
		want codes.Code
	}{
		{"not running", manager.ErrNotRunning, codes.FailedPrecondition},
		{"overloaded", manager.ErrOverloaded, codes.ResourceExhausted},
		{"invalid record", manager.ErrInvalidRecord, codes.InvalidArgument},
		{"request too large", manager.ErrRequestTooLarge, codes.InvalidArgument},
		{"unknown session", manager.ErrSessionNotFound, codes.NotFound},
		{"empty session", manager.ErrEmptySessionID, codes.InvalidArgument},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			env.runs.upsertErr = tc.err
			_, err := env.client.UpsertScalars(ctx, &corev1.UpsertScalarsRequest{
				SessionId: handle,
				Data:      []*datav1.ScalarRecord{{Key: "loss"}},
			})
			if status.Code(err) != tc.want {
				t.Fatalf("err=%v, want code %s", err, tc.want)
			}
		})
	}
}

func TestSharedRunIDKeepsSessionsSeparate(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)

	// 多 rank、resume 复用实验时会使用同一个 run_id，core 仍须各自建会话
	record := &runv1.StartRecord{Project: "smoke", Id: "shared-run"}
	rank0, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  record,
	})
	if err != nil {
		t.Fatalf("rank0 start: %v", err)
	}
	rank1, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  record,
	})
	if err != nil {
		t.Fatalf("rank1 start: %v", err)
	}
	if rank0.GetSessionId() == rank1.GetSessionId() {
		t.Fatal("same run_id must still get separate sessions")
	}
	if rank0.GetRun().GetId() != "shared-run" || rank1.GetRun().GetId() != "shared-run" {
		t.Fatalf("run_id must be echoed unchanged: %q / %q", rank0.GetRun().GetId(), rank1.GetRun().GetId())
	}

	// 一个会话结束并释放，不影响另一个
	if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{
		SessionId:    rank0.GetSessionId(),
		FinishRecord: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED},
	}); err != nil {
		t.Fatalf("finish rank0: %v", err)
	}
	if _, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{SessionId: rank0.GetSessionId()}); err != nil {
		t.Fatalf("confirm rank0: %v", err)
	}
	stats, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: rank1.GetSessionId()})
	if err != nil || stats.GetStats().GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("rank1 must stay running: resp=%v err=%v", stats, err)
	}
}

func TestRunLifecycleAndErrorMapping(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)

	if _, err := env.client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("missing start_record err=%v, want InvalidArgument", err)
	}
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("empty handle err=%v, want InvalidArgument", err)
	}
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: "missing"}); status.Code(err) != codes.NotFound {
		t.Fatalf("unknown handle err=%v, want NotFound", err)
	}

	handle := startRun(t, env, &runv1.StartRecord{})
	resp, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: handle})
	if err != nil || resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("running stats: resp=%v err=%v", resp, err)
	}
	if _, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{SessionId: handle}); status.Code(err) != codes.FailedPrecondition {
		t.Fatalf("early confirm err=%v, want FailedPrecondition", err)
	}
	if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{
		SessionId:    handle,
		FinishRecord: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED},
	}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	resp, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: handle})
	if err != nil || resp.GetStats().GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("draining stats: resp=%v err=%v, want RUNNING", resp, err)
	}
	if got := resp.GetStats().GetTotalRecords(); got != 7 {
		t.Fatalf("stats passthrough total_records = %d, want 7", got)
	}
	if _, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{SessionId: handle}); err != nil {
		t.Fatalf("ConfirmRunFinish: %v", err)
	}
	select {
	case <-env.ctrl.Done():
		t.Fatal("ConfirmRunFinish must not stop the server")
	case <-time.After(noShutdownWindow):
	}
	if _, err := env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: handle}); status.Code(err) != codes.NotFound {
		t.Fatalf("released handle err=%v, want NotFound", err)
	}
}

func TestDeliverRunFinishValidationAndFailureResponse(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)
	handle := startRun(t, env, &runv1.StartRecord{})

	if _, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{SessionId: handle}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("missing finish_record err=%v, want InvalidArgument", err)
	}

	record := &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED}
	env.runs.mu.Lock()
	env.runs.finishErr = fmt.Errorf("%w: close intake failed", manager.ErrRunFailed)
	env.runs.mu.Unlock()

	resp, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{
		SessionId:    handle,
		FinishRecord: record,
	})
	if err != nil {
		t.Fatalf("finish failure must be carried in the response, got err=%v", err)
	}
	if resp.GetSuccess() || !strings.Contains(resp.GetMessage(), "close intake") {
		t.Fatalf("finish failure response = %v, want success=false with message", resp)
	}

	env.runs.mu.Lock()
	got := env.runs.finishes[len(env.runs.finishes)-1]
	env.runs.mu.Unlock()
	if got.sessionID != handle || !proto.Equal(got.record, record) {
		t.Fatalf("finish record not forwarded: %+v", got)
	}
}

func TestConfirmRunFinishCarriesReleaseFailure(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)
	handle := startRun(t, env, &runv1.StartRecord{})
	if _, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{
		SessionId:    handle,
		FinishRecord: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED},
	}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}

	env.runs.mu.Lock()
	env.runs.confirmResult = manager.ConfirmResult{
		Message:       "Failed to finish run, but it has been saved locally.",
		ReleaseFailed: true,
	}
	env.runs.mu.Unlock()

	resp, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{SessionId: handle})
	if err != nil {
		t.Fatalf("release failure must be carried in the response, got err=%v", err)
	}
	if resp.GetSuccess() || !strings.Contains(resp.GetMessage(), "saved locally") {
		t.Fatalf("confirm response = %v, want success=false with local-copy message", resp)
	}
	if _, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: handle}); err != nil {
		t.Fatalf("session must be retained for retry: %v", err)
	}

	env.runs.mu.Lock()
	env.runs.confirmResult = manager.ConfirmResult{Message: "OK"}
	env.runs.mu.Unlock()
	resp, err = env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{SessionId: handle})
	if err != nil || !resp.GetSuccess() || resp.GetMessage() != "OK" {
		t.Fatalf("retry confirm resp=%v err=%v, want success with OK", resp, err)
	}
}

func TestConfirmRunFinishRetryableReleaseReturnsUnavailable(t *testing.T) {
	env := newTestEnv(t)
	spinup(t, env)
	ctx := callCtx(t)
	handle := startRun(t, env, &runv1.StartRecord{})
	if _, err := env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{
		SessionId:    handle,
		FinishRecord: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED},
	}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}

	env.runs.mu.Lock()
	env.runs.confirmResult = manager.ConfirmResult{
		Message:       "Failed to finish run, but it has been saved locally.",
		ReleaseFailed: true,
		Retryable:     true,
	}
	env.runs.mu.Unlock()

	_, err := env.client.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{SessionId: handle})
	if code := status.Code(err); code != codes.Unavailable {
		t.Fatalf("code = %v, want Unavailable (err=%v)", code, err)
	}
	if _, err = env.client.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{SessionId: handle}); err != nil {
		t.Fatalf("session must be retained for retry: %v", err)
	}
}

func TestManagerStatusMapsSessionCapacity(t *testing.T) {
	err := managerStatus(fmt.Errorf("start: %w", manager.ErrTooManySessions))
	if code := status.Code(err); code != codes.ResourceExhausted {
		t.Fatalf("code = %v, want ResourceExhausted", code)
	}
}

func TestManagerStatusMapsRunFailureAndTimeout(t *testing.T) {
	cases := []struct {
		name string
		err  error
		want codes.Code
	}{
		{"failed", fmt.Errorf("confirm: %w", manager.ErrRunFailed), codes.FailedPrecondition},
		{"deadline", fmt.Errorf("wait: %w", context.DeadlineExceeded), codes.DeadlineExceeded},
		{"canceled", fmt.Errorf("wait: %w", context.Canceled), codes.Canceled},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if code := status.Code(managerStatus(tc.err)); code != tc.want {
				t.Fatalf("code = %v, want %v", code, tc.want)
			}
		})
	}
}

func TestManagerStatusMapsClientErrors(t *testing.T) {
	cases := []struct {
		name string
		err  error
		want codes.Code
	}{
		{"invalid settings", manager.ErrInvalidSettings, codes.InvalidArgument},
		{"invalid start", fmt.Errorf("start: %w", manager.ErrInvalidStartRecord), codes.InvalidArgument},
		{"invalid finish", fmt.Errorf("finish: %w", manager.ErrInvalidFinishRecord), codes.InvalidArgument},
		{"missing key", client.ErrMissingAPIKey, codes.InvalidArgument},
		{"auth", fmt.Errorf("start: %w", client.ErrNotAuthenticated), codes.Unauthenticated},
		{"not found", &client.APIError{StatusCode: 404}, codes.NotFound},
		{"unavailable", &client.APIError{StatusCode: 503, Retryable: true}, codes.Unavailable},
		{"bad request", &client.APIError{StatusCode: 400}, codes.Internal},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if code := status.Code(managerStatus(tc.err)); code != tc.want {
				t.Fatalf("code = %v, want %v", code, tc.want)
			}
		})
	}
}

func TestCoreStateMapping(t *testing.T) {
	cases := []struct {
		state manager.RunState
		want  operationv1.CoreState
	}{
		{manager.RunStateRunning, operationv1.CoreState_CORE_STATE_RUNNING},
		{manager.RunStateDraining, operationv1.CoreState_CORE_STATE_RUNNING},
		{manager.RunStateDrained, operationv1.CoreState_CORE_STATE_FINISHED},
		{manager.RunStateFailed, operationv1.CoreState_CORE_STATE_FINISHED},
	}
	for _, tc := range cases {
		if got := coreState(tc.state); got != tc.want {
			t.Fatalf("coreState(%s) = %v, want %v", tc.state, got, tc.want)
		}
	}
}

type trackedTestEnv struct {
	ctrl    *server.Controller
	svc     *CoreService
	tracker *server.ConnTracker
	runs    *fakeManager
	lis     *bufconn.Listener
}

func newTrackedTestEnv(t *testing.T, grace time.Duration) *trackedTestEnv {
	t.Helper()
	tracker := server.NewConnTracker()
	g := grpc.NewServer(grpc.StatsHandler(tracker))
	ctrl := server.NewController(g, grace)
	runs := newFakeManager()
	runs.disconnectNotify = make(chan string, 16)
	svc := NewCoreService(ctrl, runs)
	tracker.SetOnDisconnect(svc.HandleConnDisconnect)
	ctrl.SetOnShutdown(svc.WaitDisconnectCleanups)
	svc.Register(g)

	lis := bufconn.Listen(bufconnSize)
	serveDone := make(chan struct{})
	go func() {
		defer close(serveDone)
		_ = g.Serve(lis)
	}()

	t.Cleanup(func() {
		ctrl.Shutdown("test-cleanup")
		<-ctrl.Done()
		_ = lis.Close()
		<-serveDone
	})
	return &trackedTestEnv{ctrl: ctrl, svc: svc, tracker: tracker, runs: runs, lis: lis}
}

func (e *trackedTestEnv) dial(t *testing.T) (*grpc.ClientConn, corev1.CoreServiceClient) {
	t.Helper()
	return dialBufconn(t, e.lis)
}

func TestConnDisconnectCleansOnlyDisconnectedClientSessions(t *testing.T) {
	env := newTrackedTestEnv(t, 2*time.Second)
	connA, clientA := env.dial(t)
	connB, clientB := env.dial(t)

	ctx := callCtx(t)
	if _, err := clientA.SpinupService(ctx, &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService A: %v", err)
	}
	startA, err := clientA.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  &runv1.StartRecord{Project: "smoke"},
	})
	if err != nil {
		t.Fatalf("DeliverRunStart A: %v", err)
	}
	startB, err := clientB.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  &runv1.StartRecord{Project: "smoke"},
	})
	if err != nil {
		t.Fatalf("DeliverRunStart B: %v", err)
	}

	// A 正常 finish + confirm 并断开连接：已释放会话不得再触发 DisconnectSession。
	if _, finErr := clientA.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{
		SessionId:    startA.GetSessionId(),
		FinishRecord: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED},
	}); finErr != nil {
		t.Fatalf("DeliverRunFinish A: %v", finErr)
	}
	if _, confErr := clientA.ConfirmRunFinish(ctx, &corev1.ConfirmRunFinishRequest{
		SessionId: startA.GetSessionId(),
	}); confErr != nil {
		t.Fatalf("ConfirmRunFinish A: %v", confErr)
	}
	if closeErr := connA.Close(); closeErr != nil {
		t.Fatalf("connA.Close: %v", closeErr)
	}

	// B 仍在运行，确认未被 A 的断连误伤。
	statsB, statsErr := clientB.GetOperationStats(ctx, &corev1.GetOperationStatsRequest{
		SessionId: startB.GetSessionId(),
	})
	if statsErr != nil || statsB.GetStats().GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("session B must stay running after A closes: resp=%v err=%v", statsB, statsErr)
	}
	if got := env.runs.disconnectedSnapshot(); len(got) != 0 {
		t.Fatalf("confirmed session A must not trigger disconnect cleanup, got %v", got)
	}

	// B 未经 finish 直接断开连接：服务端自动清理 B 的会话。
	if closeErr := connB.Close(); closeErr != nil {
		t.Fatalf("connB.Close: %v", closeErr)
	}
	select {
	case sid := <-env.runs.disconnectNotify:
		if sid != startB.GetSessionId() {
			t.Fatalf("disconnected session = %q, want %q", sid, startB.GetSessionId())
		}
	case <-time.After(2 * time.Second):
		t.Fatal("timed out waiting for disconnect cleanup of session B")
	}
}

func TestDuplicateDisconnectIsIdempotent(t *testing.T) {
	env := newTrackedTestEnv(t, 2*time.Second)
	connCtx := env.tracker.TagConn(context.Background(), &stats.ConnTagInfo{})
	env.tracker.HandleConn(connCtx, &stats.ConnBegin{})

	ctx := callCtx(t)
	if _, err := env.svc.SpinupService(ctx, &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
	if _, err := env.svc.DeliverRunStart(connCtx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  &runv1.StartRecord{Project: "smoke"},
	}); err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}

	connID, ok := server.ConnIDFromContext(connCtx)
	if !ok {
		t.Fatal("expected connID in tagged context")
	}
	env.tracker.HandleConn(connCtx, &stats.ConnEnd{})
	env.tracker.HandleConn(connCtx, &stats.ConnEnd{})
	env.svc.HandleConnDisconnect(connID)

	select {
	case <-env.runs.disconnectNotify:
	case <-time.After(2 * time.Second):
		t.Fatal("expected disconnect cleanup notification")
	}
	if err := env.svc.WaitDisconnectCleanups(ctx); err != nil {
		t.Fatalf("WaitDisconnectCleanups: %v", err)
	}
	if got := env.runs.disconnectedSnapshot(); len(got) != 1 {
		t.Fatalf("DisconnectSession called %d times, want 1", len(got))
	}
}

func TestDisconnectDuringDeliverRunStartCleansOrphanSession(t *testing.T) {
	env := newTrackedTestEnv(t, 2*time.Second)
	connCtx := env.tracker.TagConn(context.Background(), &stats.ConnTagInfo{})
	env.tracker.HandleConn(connCtx, &stats.ConnBegin{})

	ctx := callCtx(t)
	if _, err := env.svc.SpinupService(ctx, &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
	// 模拟 DeliverRunStart 执行期间底层连接突然断开的竞态。
	env.runs.onStart = func() {
		env.tracker.HandleConn(connCtx, &stats.ConnEnd{})
	}

	_, err := env.svc.DeliverRunStart(connCtx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  &runv1.StartRecord{Project: "smoke"},
	})
	if status.Code(err) != codes.Canceled {
		t.Fatalf("DeliverRunStart code = %v, want Canceled (err=%v)", status.Code(err), err)
	}

	select {
	case sid := <-env.runs.disconnectNotify:
		if sid == "" {
			t.Fatal("expected non-empty orphan session id")
		}
	case <-time.After(2 * time.Second):
		t.Fatal("timed out waiting for orphan session cleanup")
	}
}

func TestConcurrentClientsSessionLimitAndRecovery(t *testing.T) {
	env := newTrackedTestEnv(t, 2*time.Second)
	env.runs.maxSessions = 2

	conn1, client1 := env.dial(t)
	_, client2 := env.dial(t)
	_, client3 := env.dial(t)

	ctx := callCtx(t)
	if _, err := client1.SpinupService(ctx, &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
	req := func() *corev1.DeliverRunStartRequest {
		return &corev1.DeliverRunStartRequest{
			CoreSettings: validSettings(t),
			StartRecord:  &runv1.StartRecord{Project: "smoke"},
		}
	}

	if _, err := client1.DeliverRunStart(ctx, req()); err != nil {
		t.Fatalf("client1 DeliverRunStart: %v", err)
	}
	if _, err := client2.DeliverRunStart(ctx, req()); err != nil {
		t.Fatalf("client2 DeliverRunStart: %v", err)
	}
	if _, err := client3.DeliverRunStart(ctx, req()); status.Code(err) != codes.ResourceExhausted {
		t.Fatalf("client3 over-limit code = %v, want ResourceExhausted", status.Code(err))
	}

	// 客户端 1 异常断开后，会话清理释放槽位，客户端 3 可以立刻重新启动。
	if err := conn1.Close(); err != nil {
		t.Fatalf("conn1.Close: %v", err)
	}
	select {
	case <-env.runs.disconnectNotify:
	case <-time.After(2 * time.Second):
		t.Fatal("timed out waiting for client1 cleanup")
	}
	if _, err := client3.DeliverRunStart(ctx, req()); err != nil {
		t.Fatalf("client3 DeliverRunStart after disconnect cleanup: %v", err)
	}
}

func TestWaitDisconnectCleanupsCancelsSlowDrainOnShutdown(t *testing.T) {
	env := newTrackedTestEnv(t, 100*time.Millisecond)
	entered := make(chan struct{})
	canceled := make(chan struct{})
	env.runs.onDisconnect = func(ctx context.Context, _ string) {
		close(entered)
		<-ctx.Done()
		close(canceled)
	}

	conn, client := env.dial(t)
	ctx := callCtx(t)
	if _, err := client.SpinupService(ctx, &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
	if _, err := client.DeliverRunStart(ctx, &corev1.DeliverRunStartRequest{
		CoreSettings: validSettings(t),
		StartRecord:  &runv1.StartRecord{Project: "smoke"},
	}); err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	if err := conn.Close(); err != nil {
		t.Fatalf("conn.Close: %v", err)
	}

	select {
	case <-entered:
	case <-time.After(2 * time.Second):
		t.Fatal("timed out waiting for DisconnectSession to start")
	}

	env.ctrl.Shutdown("grace-timeout")
	<-env.ctrl.Done()

	select {
	case <-canceled:
	case <-time.After(2 * time.Second):
		t.Fatal("WaitDisconnectCleanups must cancel in-flight drain ctx on shutdown grace expiry")
	}
}

type notifyRunManager struct {
	*manager.Manager
	disconnected chan string
}

func (n *notifyRunManager) DisconnectSession(ctx context.Context, sessionID string) error {
	err := n.Manager.DisconnectSession(ctx, sessionID)
	n.disconnected <- sessionID
	return err
}

func TestGRPCConfirmAfterDisconnectReturnsNotFound(t *testing.T) {
	mgr := manager.New(manager.Options{})
	t.Cleanup(func() { _ = mgr.Close() })
	runs := &notifyRunManager{
		Manager:      mgr,
		disconnected: make(chan string, 4),
	}

	tracker := server.NewConnTracker()
	g := grpc.NewServer(grpc.StatsHandler(tracker))
	ctrl := server.NewController(g, 2*time.Second)
	svc := NewCoreService(ctrl, runs)
	tracker.SetOnDisconnect(svc.HandleConnDisconnect)
	ctrl.SetOnShutdown(func(ctx context.Context) error {
		_ = svc.WaitDisconnectCleanups(ctx)
		return mgr.Shutdown(ctx)
	})
	svc.Register(g)

	lis := bufconn.Listen(bufconnSize)
	serveDone := make(chan struct{})
	go func() {
		defer close(serveDone)
		_ = g.Serve(lis)
	}()
	t.Cleanup(func() {
		ctrl.Shutdown("test-cleanup")
		<-ctrl.Done()
		_ = lis.Close()
		<-serveDone
	})

	connA, clientA := dialBufconn(t, lis)
	_, clientB := dialBufconn(t, lis)

	if _, err := clientA.SpinupService(callCtx(t), &corev1.SpinupServiceRequest{}); err != nil {
		t.Fatalf("SpinupService: %v", err)
	}
	settings := validSettings(t)
	settings.Mode = settingsv1.CoreMode_CORE_MODE_LOCAL
	startResp, err := clientA.DeliverRunStart(callCtx(t), &corev1.DeliverRunStartRequest{
		CoreSettings: settings,
		StartRecord:  &runv1.StartRecord{Project: "smoke"},
	})
	if err != nil {
		t.Fatalf("DeliverRunStart: %v", err)
	}
	sessionID := startResp.GetSessionId()

	// 先完成 DeliverRunFinish 但尚未 Confirm 时连接断开，触发断连排空清理。
	if _, err = clientA.DeliverRunFinish(callCtx(t), &corev1.DeliverRunFinishRequest{
		SessionId:    sessionID,
		FinishRecord: &runv1.FinishRecord{State: runv1.RunState_RUN_STATE_FINISHED},
	}); err != nil {
		t.Fatalf("DeliverRunFinish: %v", err)
	}
	if err = connA.Close(); err != nil {
		t.Fatalf("connA.Close: %v", err)
	}

	select {
	case sid := <-runs.disconnected:
		if sid != sessionID {
			t.Fatalf("disconnected session = %q, want %q", sid, sessionID)
		}
	case <-time.After(callTimeout):
		t.Fatal("timed out waiting for DisconnectSession after connA close")
	}
	if got := mgr.SessionCount(); got != 0 {
		t.Fatalf("SessionCount after disconnect cleanup = %d, want 0", got)
	}

	// 重连后重发 ConfirmRunFinish 与 UpsertScalars 都必须明确失败（codes.NotFound），不得假成功。
	if _, err = clientB.ConfirmRunFinish(callCtx(t), &corev1.ConfirmRunFinishRequest{SessionId: sessionID}); status.Code(err) != codes.NotFound {
		t.Fatalf("ConfirmRunFinish after disconnect code = %v, want NotFound (err=%v)", status.Code(err), err)
	}
	if _, err = clientB.UpsertScalars(callCtx(t), &corev1.UpsertScalarsRequest{
		SessionId: sessionID,
		Data:      []*datav1.ScalarRecord{{Key: "loss"}},
	}); status.Code(err) != codes.NotFound {
		t.Fatalf("UpsertScalars after disconnect code = %v, want NotFound (err=%v)", status.Code(err), err)
	}
}
