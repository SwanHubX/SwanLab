package service

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"net"
	"net/http"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"

	"google.golang.org/grpc"
	"google.golang.org/grpc/codes"
	"google.golang.org/grpc/credentials/insecure"
	"google.golang.org/grpc/status"
	"google.golang.org/grpc/test/bufconn"
	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/emptypb"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/server"
	"github.com/swanhubx/swanlab/core/internal/store"
	"github.com/swanhubx/swanlab/core/internal/syncer"
	corev1 "github.com/swanhubx/swanlab/core/proto/swanlab/grpc/core/v1"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

type fakeSyncManager struct {
	mu           sync.Mutex
	next         int
	sessions     map[string]string
	starts       []syncer.StartOptions
	flushed      []string
	confirmed    []string
	disconnected []string

	startErr   error
	flushErr   error
	statsErr   error
	confirmErr error

	flushResult   syncer.FlushResult
	confirmResult syncer.ConfirmResult

	onStart          func()
	onDisconnect     func(context.Context, string)
	disconnectNotify chan string
}

func newFakeSyncManager() *fakeSyncManager {
	return &fakeSyncManager{
		sessions:      make(map[string]string),
		flushResult:   syncer.FlushResult{Path: "/ws/proj/run-1"},
		confirmResult: syncer.ConfirmResult{Success: true, Message: "success"},
	}
}

func (f *fakeSyncManager) DeliverSyncStart(
	_ context.Context,
	opts syncer.StartOptions,
) (string, error) {
	f.mu.Lock()
	hook := f.onStart
	f.mu.Unlock()
	if hook != nil {
		hook()
	}

	f.mu.Lock()
	defer f.mu.Unlock()
	if f.startErr != nil {
		return "", f.startErr
	}
	f.next++
	id := fmt.Sprintf("sync-%d", f.next)
	f.starts = append(f.starts, opts)
	f.sessions[id] = "started"
	return id, nil
}

func (f *fakeSyncManager) DeliverSyncFlush(_ context.Context, sessionID string) (syncer.FlushResult, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.flushErr != nil {
		return syncer.FlushResult{}, f.flushErr
	}
	if _, ok := f.sessions[sessionID]; !ok {
		return syncer.FlushResult{}, syncer.ErrSessionNotFound
	}
	f.sessions[sessionID] = "running"
	f.flushed = append(f.flushed, sessionID)
	return f.flushResult, nil
}

func (f *fakeSyncManager) Stats(sessionID string) (*operationv1.OperationStats, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.statsErr != nil {
		return nil, f.statsErr
	}
	if _, ok := f.sessions[sessionID]; !ok {
		return nil, syncer.ErrSessionNotFound
	}
	return &operationv1.OperationStats{
		State:           operationv1.CoreState_CORE_STATE_FINISHED,
		TotalRecords:    4,
		UploadedRecords: 4,
	}, nil
}

func (f *fakeSyncManager) ConfirmSyncFinish(_ context.Context, sessionID string) (syncer.ConfirmResult, error) {
	f.mu.Lock()
	defer f.mu.Unlock()
	if f.confirmErr != nil {
		return syncer.ConfirmResult{}, f.confirmErr
	}
	if _, ok := f.sessions[sessionID]; !ok {
		return syncer.ConfirmResult{}, syncer.ErrSessionNotFound
	}
	f.confirmed = append(f.confirmed, sessionID)
	delete(f.sessions, sessionID)
	return f.confirmResult, nil
}

func (f *fakeSyncManager) Disconnect(ctx context.Context, sessionID string) error {
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
	f.mu.Unlock()

	if notify != nil {
		notify <- sessionID
	}
	return nil
}

type syncTestEnv struct {
	syncClient corev1.CoreSyncServiceClient
	coreClient corev1.CoreServiceClient
	ctrl       *server.Controller
	svc        *SyncService
	mgr        *fakeSyncManager
	lis        *bufconn.Listener
}

func dialSyncBufconn(
	t *testing.T,
	lis *bufconn.Listener,
) (*grpc.ClientConn, corev1.CoreServiceClient, corev1.CoreSyncServiceClient) {
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
	return conn, corev1.NewCoreServiceClient(conn), corev1.NewCoreSyncServiceClient(conn)
}

func newSyncTestEnv(t *testing.T, tracked bool) *syncTestEnv {
	t.Helper()
	var opts []grpc.ServerOption
	var tracker *server.ConnTracker
	if tracked {
		tracker = server.NewConnTracker()
		opts = append(opts, grpc.StatsHandler(tracker))
	}
	g := grpc.NewServer(opts...)
	ctrl := server.NewController(g, 2*time.Second)
	mgr := newFakeSyncManager()
	svc := NewSyncService(ctrl, mgr)
	svc.Register(g)
	if tracker != nil {
		tracker.SetOnDisconnect(svc.HandleConnDisconnect)
	}
	lis := bufconn.Listen(bufconnSize)
	go func() { _ = g.Serve(lis) }()
	t.Cleanup(func() {
		_ = svc.WaitDisconnectCleanups(context.Background())
		ctrl.Shutdown("test-cleanup")
		<-ctrl.Done()
	})
	_, coreCli, syncCli := dialSyncBufconn(t, lis)
	return &syncTestEnv{
		syncClient: syncCli,
		coreClient: coreCli,
		ctrl:       ctrl,
		svc:        svc,
		mgr:        mgr,
		lis:        lis,
	}
}

func TestSyncServiceFullLifecycleWithFake(t *testing.T) {
	env := newSyncTestEnv(t, false)
	ctx := callCtx(t)

	startResp, err := env.syncClient.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
		CoreSettings: validSettings(t),
		Workspace:    "ws-override",
		Project:      "proj-override",
		Id:           "run-override",
	})
	if err != nil || !startResp.GetSuccess() {
		t.Fatalf("DeliverSyncStart: resp=%v err=%v", startResp, err)
	}

	flushResp, err := env.syncClient.DeliverSyncFlush(ctx, &emptypb.Empty{})
	if err != nil || !flushResp.GetSuccess() || flushResp.GetPath() != "/ws/proj/run-1" {
		t.Fatalf("DeliverSyncFlush: resp=%v err=%v", flushResp, err)
	}

	statsResp, err := env.syncClient.GetOperationStats(ctx, &emptypb.Empty{})
	if err != nil || !statsResp.GetSuccess() || statsResp.GetStats().GetTotalRecords() != 4 {
		t.Fatalf("GetOperationStats: resp=%v err=%v", statsResp, err)
	}

	confirmResp, err := env.syncClient.ConfirmSyncFinish(ctx, &emptypb.Empty{})
	if err != nil || !confirmResp.GetSuccess() {
		t.Fatalf("ConfirmSyncFinish: resp=%v err=%v", confirmResp, err)
	}

	env.mgr.mu.Lock()
	if len(env.mgr.starts) != 1 || env.mgr.starts[0].Workspace != "ws-override" || env.mgr.starts[0].RunID != "run-override" {
		env.mgr.mu.Unlock()
		t.Fatalf("unexpected starts recorded: %+v", env.mgr.starts)
	}
	env.mgr.mu.Unlock()

	// ConfirmSyncFinish 后会话已解绑，再次查询统计返回 Success=false，Flush 返回 FailedPrecondition。
	postStats, err := env.syncClient.GetOperationStats(ctx, &emptypb.Empty{})
	if err != nil || postStats.GetSuccess() {
		t.Fatalf("post-confirm GetOperationStats = (%v, %v), want Success=false", postStats, err)
	}
	if _, err := env.syncClient.DeliverSyncFlush(ctx, &emptypb.Empty{}); status.Code(err) != codes.FailedPrecondition {
		t.Fatalf("post-confirm DeliverSyncFlush code = %v, want FailedPrecondition", status.Code(err))
	}
}

func TestSyncServiceProtocolPreconditions(t *testing.T) {
	env := newSyncTestEnv(t, false)
	ctx := callCtx(t)

	if _, err := env.syncClient.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{}); status.Code(err) != codes.InvalidArgument {
		t.Fatalf("nil CoreSettings code = %v, want InvalidArgument", status.Code(err))
	}
	if _, err := env.syncClient.DeliverSyncFlush(ctx, &emptypb.Empty{}); status.Code(err) != codes.FailedPrecondition {
		t.Fatalf("flush before start code = %v, want FailedPrecondition", status.Code(err))
	}
	if _, err := env.syncClient.ConfirmSyncFinish(ctx, &emptypb.Empty{}); status.Code(err) != codes.FailedPrecondition {
		t.Fatalf("confirm before start code = %v, want FailedPrecondition", status.Code(err))
	}
	statsResp, err := env.syncClient.GetOperationStats(ctx, &emptypb.Empty{})
	if err != nil || statsResp.GetSuccess() {
		t.Fatalf("stats before start = (%v, %v), want Success=false", statsResp, err)
	}

	if resp, startErr := env.syncClient.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
		CoreSettings: validSettings(t),
	}); startErr != nil || !resp.GetSuccess() {
		t.Fatalf("DeliverSyncStart: resp=%v err=%v", resp, startErr)
	}
	if _, dupErr := env.syncClient.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
		CoreSettings: validSettings(t),
	}); status.Code(dupErr) != codes.FailedPrecondition {
		t.Fatalf("duplicate DeliverSyncStart code = %v, want FailedPrecondition", status.Code(dupErr))
	}

	env.ctrl.Lifecycle().BeginStopping()
	if _, stopErr := env.syncClient.DeliverSyncFlush(ctx, &emptypb.Empty{}); status.Code(stopErr) != codes.FailedPrecondition {
		t.Fatalf("flush while stopping code = %v, want FailedPrecondition", status.Code(stopErr))
	}
}

func TestSyncServiceErrorMappings(t *testing.T) {
	t.Run("start domain rejection vs too many sessions", func(t *testing.T) {
		env := newSyncTestEnv(t, false)
		ctx := callCtx(t)

		for _, domainErr := range []error{
			fmt.Errorf("%w: bad dir", syncer.ErrInvalidSettings),
			fmt.Errorf("%w: missing start", syncer.ErrInvalidStartRecord),
			fmt.Errorf("%w: reversed step", syncer.ErrInvalidSyncRange),
		} {
			env.mgr.startErr = domainErr
			resp, err := env.syncClient.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
				CoreSettings: validSettings(t),
			})
			if err != nil {
				t.Fatalf("DeliverSyncStart unexpected gRPC error: %v", err)
			}
			if resp.GetSuccess() || !strings.Contains(resp.GetMessage(), domainErr.Error()) {
				t.Fatalf("DeliverSyncStart = %+v, want Success=false containing %q", resp, domainErr.Error())
			}
		}

		env.mgr.startErr = syncer.ErrTooManySessions
		_, err := env.syncClient.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
			CoreSettings: validSettings(t),
		})
		if status.Code(err) != codes.ResourceExhausted {
			t.Fatalf("ErrTooManySessions code = %v, want ResourceExhausted", status.Code(err))
		}
	})

	t.Run("flush failure vs confirm soft failure unbinds session", func(t *testing.T) {
		env := newSyncTestEnv(t, false)
		ctx := callCtx(t)

		if _, err := env.syncClient.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
			CoreSettings: validSettings(t),
		}); err != nil {
			t.Fatalf("DeliverSyncStart: %v", err)
		}

		env.mgr.flushErr = errors.New("403 forbidden")
		flushResp, err := env.syncClient.DeliverSyncFlush(ctx, &emptypb.Empty{})
		if err != nil {
			t.Fatalf("DeliverSyncFlush unexpected gRPC error: %v", err)
		}
		if flushResp.GetSuccess() || !strings.Contains(flushResp.GetMessage(), "403 forbidden") {
			t.Fatalf("DeliverSyncFlush = %+v, want Success=false", flushResp)
		}

		env.mgr.flushErr = syncer.ErrInvalidState
		if _, stErr := env.syncClient.DeliverSyncFlush(ctx, &emptypb.Empty{}); status.Code(stErr) != codes.FailedPrecondition {
			t.Fatalf("flush ErrInvalidState code = %v, want FailedPrecondition", status.Code(stErr))
		}

		env.mgr.flushErr = nil
		env.mgr.confirmResult = syncer.ConfirmResult{Success: false, Message: "stop failed"}
		confirmResp, err := env.syncClient.ConfirmSyncFinish(ctx, &emptypb.Empty{})
		if err != nil {
			t.Fatalf("ConfirmSyncFinish unexpected gRPC error: %v", err)
		}
		if confirmResp.GetSuccess() || confirmResp.GetMessage() != "stop failed" {
			t.Fatalf("ConfirmSyncFinish = %+v, want Success=false stop failed", confirmResp)
		}
		if _, afterErr := env.syncClient.ConfirmSyncFinish(ctx, &emptypb.Empty{}); status.Code(afterErr) != codes.FailedPrecondition {
			t.Fatalf("ConfirmSyncFinish after unbind code = %v, want FailedPrecondition", status.Code(afterErr))
		}
	})
}

func TestSyncServiceDisconnectIsolatesClientsAndCancelsPendingStart(t *testing.T) {
	t.Run("disconnect cleans only disconnected conn", func(t *testing.T) {
		env := newSyncTestEnv(t, true)
		ctx := callCtx(t)

		disconnectedCh := make(chan string, 4)
		env.mgr.mu.Lock()
		env.mgr.disconnectNotify = disconnectedCh
		env.mgr.mu.Unlock()

		connA, _, syncA := dialSyncBufconn(t, env.lis)
		if _, err := syncA.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
			CoreSettings: validSettings(t),
		}); err != nil {
			t.Fatalf("A DeliverSyncStart: %v", err)
		}

		_, _, syncB := dialSyncBufconn(t, env.lis)
		if _, err := syncB.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
			CoreSettings: validSettings(t),
		}); err != nil {
			t.Fatalf("B DeliverSyncStart: %v", err)
		}

		if err := connA.Close(); err != nil {
			t.Fatalf("close connA: %v", err)
		}
		select {
		case got := <-disconnectedCh:
			if got != "sync-1" {
				t.Fatalf("disconnected session = %q, want sync-1", got)
			}
		case <-time.After(callTimeout):
			t.Fatal("timed out waiting for connA sync disconnect cleanup")
		}

		if _, err := syncB.DeliverSyncFlush(ctx, &emptypb.Empty{}); err != nil {
			t.Fatalf("B DeliverSyncFlush after A disconnect: %v", err)
		}
		if resp, err := syncB.ConfirmSyncFinish(ctx, &emptypb.Empty{}); err != nil || !resp.GetSuccess() {
			t.Fatalf("B ConfirmSyncFinish = (%v, %v), want Success=true", resp, err)
		}
	})

	t.Run("disconnect during DeliverSyncStart cleans orphan session", func(t *testing.T) {
		env := newSyncTestEnv(t, true)
		conn, _, syncCli := dialSyncBufconn(t, env.lis)

		enteredStart := make(chan struct{})
		releaseStart := make(chan struct{})
		disconnectedCh := make(chan string, 2)

		env.mgr.mu.Lock()
		env.mgr.disconnectNotify = disconnectedCh
		env.mgr.onStart = func() {
			close(enteredStart)
			<-releaseStart
		}
		env.mgr.mu.Unlock()

		errCh := make(chan error, 1)
		go func() {
			_, err := syncCli.DeliverSyncStart(context.Background(), &corev1.DeliverSyncStartRequest{
				CoreSettings: validSettings(t),
			})
			errCh <- err
		}()

		<-enteredStart
		if err := conn.Close(); err != nil {
			t.Fatalf("conn.Close: %v", err)
		}
		close(releaseStart)

		if err := <-errCh; status.Code(err) != codes.Canceled {
			t.Fatalf("DeliverSyncStart after disconnect code = %v, want Canceled (err=%v)", status.Code(err), err)
		}
		select {
		case sid := <-disconnectedCh:
			if sid == "" {
				t.Fatal("expected non-empty cleaned session ID")
			}
		case <-time.After(callTimeout):
			t.Fatal("timed out waiting for orphan session disconnect")
		}
	})

	t.Run("WaitDisconnectCleanups cancels slow cleanup on timeout", func(t *testing.T) {
		env := newSyncTestEnv(t, true)
		ctx := callCtx(t)

		cleanupStarted := make(chan struct{})
		cleanupCanceled := make(chan struct{})
		env.mgr.mu.Lock()
		env.mgr.onDisconnect = func(dctx context.Context, _ string) {
			close(cleanupStarted)
			<-dctx.Done()
			close(cleanupCanceled)
		}
		env.mgr.mu.Unlock()

		conn, _, syncCli := dialSyncBufconn(t, env.lis)
		if _, err := syncCli.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
			CoreSettings: validSettings(t),
		}); err != nil {
			t.Fatalf("DeliverSyncStart: %v", err)
		}
		_ = conn.Close()
		<-cleanupStarted

		waitCtx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
		defer cancel()
		if err := env.svc.WaitDisconnectCleanups(waitCtx); !errors.Is(err, context.DeadlineExceeded) {
			t.Fatalf("WaitDisconnectCleanups = %v, want DeadlineExceeded", err)
		}
		<-cleanupCanceled
	})
}

func TestSyncServiceEndToEndWithRealSyncer(t *testing.T) {
	cloud := clienttest.NewCloud(t)

	var (
		mu         sync.Mutex
		metricReqs int
		stopState  string
	)
	cloud.ServeJSON(http.MethodGet, "/api/auth/verify", http.StatusOK, map[string]any{
		"uid":       1,
		"username":  "alice",
		"createdAt": "2026-01-01T00:00:00Z",
	})
	cloud.ServeJSON(http.MethodGet, "/api/project/ws1/proj-sync", http.StatusOK, map[string]any{
		"cuid":       "proj-1",
		"name":       "proj-sync",
		"path":       "/@alice/proj-sync",
		"visibility": "PRIVATE",
		"version":    1,
		"_count":     map[string]any{"experiments": 1},
		"group":      map[string]any{"username": "alice"},
	})
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/proj-sync/experiment", http.StatusCreated, map[string]any{
		"cuid":      "exp-1",
		"slug":      "exp-slug",
		"name":      "exp1",
		"createdAt": "2026-01-01T00:00:00Z",
	})
	cloud.ServeJSON(http.MethodPost, "/api/projects/alice/proj-sync/series", http.StatusOK, map[string]any{"ok": true})
	cloud.Handle(http.MethodPut, "/api/project/alice/proj-sync/runs/exp-1/state", func(w http.ResponseWriter, r *http.Request) {
		var body map[string]any
		if err := json.NewDecoder(r.Body).Decode(&body); err == nil {
			if st, ok := body["state"].(string); ok {
				mu.Lock()
				stopState = st
				mu.Unlock()
			}
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, _ *http.Request) {
		mu.Lock()
		metricReqs++
		mu.Unlock()
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
	})

	runDir := t.TempDir()
	runFile := filepath.Join(runDir, "run-source.swanlab")
	w, err := store.Create(runFile)
	if err != nil {
		t.Fatalf("store.Create: %v", err)
	}
	ts := timestamppb.Now()
	recs := []*recordv1.Record{
		{Num: 1, Timestamp: ts, RecordType: &recordv1.Record_Start{Start: &runv1.StartRecord{
			Id: "orig-run", Workspace: "ws1", Project: "proj1", Name: "exp1",
		}}},
		{Num: 2, Timestamp: ts, RecordType: &recordv1.Record_Column{Column: &columnv1.ColumnRecord{
			ColumnKey: "loss", ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR, ColumnClass: columnv1.ColumnClass_COLUMN_CLASS_CUSTOM,
		}}},
		{Num: 3, Timestamp: ts, RecordType: &recordv1.Record_Scalar{Scalar: &datav1.ScalarRecord{
			Key: "loss", Step: 1, Type: columnv1.ColumnType_COLUMN_TYPE_SCALAR, Value: &datav1.ScalarValue{Number: 0.25}, Timestamp: ts,
		}}},
		{Num: 4, Timestamp: ts, RecordType: &recordv1.Record_Log{Log: &terminalv1.LogRecord{
			Level: terminalv1.LogLevel_LOG_LEVEL_INFO, Epoch: 1, Timestamp: ts, Line: "hello sync service",
		}}},
		{Num: 5, Timestamp: ts, RecordType: &recordv1.Record_Finish{Finish: &runv1.FinishRecord{
			State: runv1.RunState_RUN_STATE_FINISHED,
		}}},
	}
	for _, rec := range recs {
		raw, mErr := proto.Marshal(rec)
		if mErr != nil {
			t.Fatalf("proto.Marshal: %v", mErr)
		}
		if wErr := w.Write(raw); wErr != nil {
			t.Fatalf("store.Write: %v", wErr)
		}
	}
	if cErr := w.Close(); cErr != nil {
		t.Fatalf("writer.Close: %v", cErr)
	}

	tracker := server.NewConnTracker()
	g := grpc.NewServer(grpc.StatsHandler(tracker))
	ctrl := server.NewController(g, 2*time.Second)
	coreMgr := newFakeManager()
	coreSvc := NewCoreService(ctrl, coreMgr)
	syncMgr := syncer.NewManager(syncer.Options{HTTP: client.Options{Version: "test"}})
	syncSvc := NewSyncService(ctrl, syncMgr)
	tracker.SetOnDisconnect(func(connID uint64) {
		coreSvc.HandleConnDisconnect(connID)
		syncSvc.HandleConnDisconnect(connID)
	})
	coreSvc.Register(g)
	syncSvc.Register(g)

	lis := bufconn.Listen(bufconnSize)
	go func() { _ = g.Serve(lis) }()
	t.Cleanup(func() {
		_ = syncSvc.WaitDisconnectCleanups(context.Background())
		_ = syncMgr.Shutdown(context.Background())
		ctrl.Shutdown("e2e-cleanup")
		<-ctrl.Done()
	})

	_, coreCli, syncCli := dialSyncBufconn(t, lis)
	ctx := callCtx(t)
	if _, spErr := coreCli.SpinupService(ctx, &corev1.SpinupServiceRequest{}); spErr != nil {
		t.Fatalf("SpinupService: %v", spErr)
	}

	settings := validSettings(t)
	settings.RunDir = runDir
	settings.ApiHost = cloud.URL()
	settings.ApiKey = "test-api-key"

	startResp, err := syncCli.DeliverSyncStart(ctx, &corev1.DeliverSyncStartRequest{
		CoreSettings: settings,
		Project:      "proj-sync",
		Id:           "target-run",
	})
	if err != nil || !startResp.GetSuccess() {
		t.Fatalf("DeliverSyncStart: resp=%v err=%v", startResp, err)
	}

	flushResp, err := syncCli.DeliverSyncFlush(ctx, &emptypb.Empty{})
	if err != nil || !flushResp.GetSuccess() || flushResp.GetPath() != "/@alice/proj-sync/exp-slug" {
		t.Fatalf("DeliverSyncFlush: resp=%v err=%v", flushResp, err)
	}

	statsResp, err := syncCli.GetOperationStats(ctx, &emptypb.Empty{})
	if err != nil || !statsResp.GetSuccess() || statsResp.GetStats() == nil {
		t.Fatalf("GetOperationStats: resp=%v err=%v", statsResp, err)
	}

	confirmResp, err := syncCli.ConfirmSyncFinish(ctx, &emptypb.Empty{})
	if err != nil || !confirmResp.GetSuccess() {
		t.Fatalf("ConfirmSyncFinish: resp=%v err=%v", confirmResp, err)
	}

	mu.Lock()
	defer mu.Unlock()
	if stopState != "FINISHED" {
		t.Fatalf("stopState = %q, want FINISHED", stopState)
	}
	if metricReqs == 0 {
		t.Fatal("expected uploaded metrics/log batches on mock cloud")
	}
}
