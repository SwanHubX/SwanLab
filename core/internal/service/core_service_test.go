package service

import (
	"context"
	"fmt"
	"net"
	"sync"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials/insecure"
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

// fakeManager 只实现 service 用到的编排接口，用于验证 RPC 契约与错误映射。
type fakeManager struct {
	mu       sync.Mutex
	next     int
	sessions map[string]manager.RunState
	starts   []fakeStart

	startErr   error
	finishErr  error
	stateErr   error
	confirmErr error
	upsertErr  error

	upserts []fakeUpsert
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
	return &fakeManager{sessions: make(map[string]manager.RunState)}
}

func (f *fakeManager) DeliverRunStart(
	_ context.Context,
	settings config.Settings,
	record *runv1.StartRecord,
) (*manager.StartResult, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	f.starts = append(f.starts, fakeStart{settings: settings, record: record})
	if f.startErr != nil {
		return nil, f.startErr
	}
	f.next++
	sessionID := fmt.Sprintf("session-%d", f.next)
	f.sessions[sessionID] = manager.RunStateRunning
	return &manager.StartResult{
		SessionID:     sessionID,
		Run:           record,
		Path:          "/alice/demo/run-1",
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

func (f *fakeManager) UpsertLogs(_ context.Context, sessionID string, records []*terminalv1.LogRecord) error {
	return f.recordUpsert("logs", sessionID, len(records))
}

func (f *fakeManager) UpsertSaves(_ context.Context, sessionID string, records []*savev1.SaveRecord) error {
	return f.recordUpsert("saves", sessionID, len(records))
}

func (f *fakeManager) Finish(sessionID string) error {
	f.mu.Lock()
	defer f.mu.Unlock()
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

func (f *fakeManager) Confirm(_ context.Context, sessionID string) error {
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.confirmErr != nil {
		return f.confirmErr
	}
	if sessionID == "" {
		return manager.ErrEmptySessionID
	}
	state, ok := f.sessions[sessionID]
	if !ok {
		return manager.ErrSessionNotFound
	}
	switch state {
	case manager.RunStateRunning:
		return manager.ErrRunNotFinished
	case manager.RunStateFailed:
		return manager.ErrRunFailed
	}
	delete(f.sessions, sessionID)
	return nil
}

type testEnv struct {
	client corev1.CoreServiceClient
	ctrl   *server.Controller
	runs   *fakeManager
}

func newTestEnv(t *testing.T) *testEnv {
	t.Helper()
	g := grpc.NewServer()
	ctrl := server.NewController(g, 2*time.Second)
	runs := newFakeManager()
	NewCoreService(ctrl, runs).Register(g)
	lis := bufconn.Listen(bufconnSize)
	go func() { _ = g.Serve(lis) }()
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
	t.Cleanup(func() {
		_ = conn.Close()
		ctrl.Shutdown("test-cleanup")
		<-ctrl.Done()
	})
	return &testEnv{client: corev1.NewCoreServiceClient(conn), ctrl: ctrl, runs: runs}
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
	}
	got, err := settingsFromProto(in)
	if err != nil {
		t.Fatalf("settingsFromProto: %v", err)
	}
	want := config.Settings{
		RunID:          "run-1",
		RunDir:         "/tmp/run",
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
	}
	if got != want {
		t.Fatalf("settings = %+v, want %+v", got, want)
	}

	if _, err := settingsFromProto(nil); err == nil {
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
			if _, err := settingsFromProto(clone); err == nil {
				t.Fatal("expected an error")
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
	if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{SessionId: rank0.GetSessionId()}); err != nil {
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
	if _, err = env.client.DeliverRunFinish(ctx, &corev1.DeliverRunFinishRequest{SessionId: handle}); err != nil {
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
