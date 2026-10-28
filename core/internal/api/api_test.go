package api

import (
	"context"
	"net/http"
	"strings"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
)

// TestCreateExperimentCreatesNew 覆盖 201 分支：新建实验以「发起请求的当前 UTC 时间」
// 作为 createdAt，且缺失字段从请求体中剔除。
func TestCreateExperimentCreatesNew(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusCreated, map[string]any{
		"cuid": "exp-1", "name": "demo", "slug": "exp-1",
	})

	svc, clock := newTestClient(t, cloud, client.Options{})
	experiment, isNew, err := svc.CreateOrResumeExperiment(context.Background(), CreateExperimentParams{
		Username: "alice",
		Project:  "demo",
		Name:     "demo",
		Color:    "#FF0000",
		Tags:     []string{"tag-a", "tag-b"},
		Resume:   ResumeAllow,
	})
	if err != nil {
		t.Fatalf("create experiment: %v", err)
	}
	if !isNew {
		t.Fatal("201 must be reported as a new experiment")
	}
	wantCreatedAt := "2026-03-10T19:02:03.123400Z"
	if experiment.CreatedAt != wantCreatedAt {
		t.Fatalf("createdAt = %q, want %q", experiment.CreatedAt, wantCreatedAt)
	}

	reqs := cloud.RequestsFor(http.MethodPost, "/api/project/alice/demo/experiment")
	if len(reqs) != 1 {
		t.Fatalf("create requests = %d, want 1", len(reqs))
	}
	body := decodeJSONBody(t, reqs[0].Body)
	if body["createdAt"] != wantCreatedAt {
		t.Fatalf("body createdAt = %v, want %q", body["createdAt"], wantCreatedAt)
	}
	colors, _ := body["colors"].([]any)
	if len(colors) != 2 || colors[0] != "#FF0000" || colors[1] != "#FF0000" {
		t.Fatalf("body colors = %v, want duplicated color", body["colors"])
	}
	labels, _ := body["labels"].([]any)
	if len(labels) != 2 {
		t.Fatalf("body labels = %v, want 2 entries", body["labels"])
	}
	for _, absent := range []string{"description", "job", "cluster", "cuid"} {
		if _, ok := body[absent]; ok {
			t.Fatalf("body must omit empty %q, got %v", absent, body[absent])
		}
	}
	if clock.Now().Sub(clienttest.FixedNow) != 0 {
		t.Fatalf("createdAt must be taken at request time, clock advanced to %v", clock.Now())
	}
	if got := cloud.CountFor(http.MethodGet, "/api/project/alice/demo/runs/exp-1"); got != 0 {
		t.Fatalf("new experiment must not query createdAt again, got %d requests", got)
	}
}

// TestCreateExperimentExistingUsesResponseCreatedAt 覆盖 200 分支：后端返回的
// createdAt 优先，不再额外查询。
func TestCreateExperimentExistingUsesResponseCreatedAt(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusOK, map[string]any{
		"cuid": "exp-1", "createdAt": "2026-01-02T03:04:05.000000Z",
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	experiment, isNew, err := svc.CreateOrResumeExperiment(context.Background(), CreateExperimentParams{
		Username: "alice",
		Project:  "demo",
		Name:     "demo",
		RunID:    "run-1",
		Resume:   ResumeAllow,
	})
	if err != nil {
		t.Fatalf("create experiment: %v", err)
	}
	if isNew {
		t.Fatal("200 must be reported as an existing experiment")
	}
	if experiment.CreatedAt != "2026-01-02T03:04:05.000000Z" {
		t.Fatalf("createdAt = %q, want backend value", experiment.CreatedAt)
	}
	body := decodeJSONBody(t, cloud.RequestsFor(http.MethodPost, "/api/project/alice/demo/experiment")[0].Body)
	if body["cuid"] != "run-1" {
		t.Fatalf("body cuid = %v, want run-1", body["cuid"])
	}
}

// TestCreateExperimentFallsBackToRunQuery 覆盖旧后端：POST 未回传 createdAt 时回退查询。
func TestCreateExperimentFallsBackToRunQuery(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusOK, map[string]any{"cuid": "exp-1"})
	cloud.ServeJSON(http.MethodGet, "/api/project/alice/demo/runs/exp-1", http.StatusOK, map[string]any{
		"cuid": "exp-1", "createdAt": "2025-12-31T00:00:00.000Z",
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	experiment, _, err := svc.CreateOrResumeExperiment(context.Background(), CreateExperimentParams{
		Username: "alice", Project: "demo", Name: "demo",
	})
	if err != nil {
		t.Fatalf("create experiment: %v", err)
	}
	if experiment.CreatedAt != "2025-12-31T00:00:00.000Z" {
		t.Fatalf("createdAt = %q, want fallback value", experiment.CreatedAt)
	}
	if got := cloud.CountFor(http.MethodGet, "/api/project/alice/demo/runs/exp-1"); got != 1 {
		t.Fatalf("fallback queries = %d, want 1", got)
	}
}

// TestCreateExperimentFailsWithoutCreatedAt 覆盖「旧旧后端」：两条路径都拿不到
// createdAt 时明确报错并提示升级，而不是静默用当前时间兜底。
func TestCreateExperimentFailsWithoutCreatedAt(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, "/api/project/alice/demo/experiment", http.StatusOK, map[string]any{"cuid": "exp-1"})
	cloud.ServeJSON(http.MethodGet, "/api/project/alice/demo/runs/exp-1", http.StatusOK, map[string]any{"cuid": "exp-1"})

	svc, _ := newTestClient(t, cloud, client.Options{})
	_, _, err := svc.CreateOrResumeExperiment(context.Background(), CreateExperimentParams{
		Username: "alice", Project: "demo", Name: "demo",
	})
	if err == nil {
		t.Fatal("expected error when backend never returns createdAt")
	}
	if !strings.Contains(err.Error(), "upgrade swanlab-server") {
		t.Fatalf("error = %v, want upgrade hint", err)
	}
}

// TestCreateExperimentResumeMustRequiresExisting 覆盖 resume=must：
// 实验不存在时报错且不发起创建请求。
func TestCreateExperimentResumeMustRequiresExisting(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, "/api/project/alice/demo/runs/missing", http.StatusNotFound, map[string]string{
		"code": "not_found", "message": "run not found",
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	_, _, err := svc.CreateOrResumeExperiment(context.Background(), CreateExperimentParams{
		Username: "alice", Project: "demo", Name: "demo", RunID: "missing", Resume: ResumeMust,
	})
	if err == nil {
		t.Fatal("expected error for missing experiment with resume=must")
	}
	if !strings.Contains(err.Error(), "does not exist in project demo") {
		t.Fatalf("error = %v, want explicit missing-experiment message", err)
	}
	if got := cloud.CountFor(http.MethodPost, "/api/project/alice/demo/experiment"); got != 0 {
		t.Fatalf("create requests = %d, want 0", got)
	}
}

// TestCreateExperimentResumeMustRequiresRunID 覆盖 resume=must 时缺少 run_id。
func TestCreateExperimentResumeMustRequiresRunID(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	svc, _ := newTestClient(t, cloud, client.Options{})
	if _, _, err := svc.CreateOrResumeExperiment(context.Background(), CreateExperimentParams{
		Username: "alice", Project: "demo", Name: "demo", Resume: ResumeMust,
	}); err == nil {
		t.Fatal("expected error when run_id is missing with resume=must")
	}
}

// TestGetOrCreateProjectExisting 覆盖项目已存在：只查询不创建。
func TestGetOrCreateProjectExisting(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, "/api/project/alice/demo", http.StatusOK, map[string]any{
		"cuid": "proj-1", "name": "demo", "visibility": "PRIVATE",
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	project, err := svc.GetOrCreateProject(context.Background(), "alice", "demo", false)
	if err != nil {
		t.Fatalf("get project: %v", err)
	}
	if project.CUID != "proj-1" {
		t.Fatalf("project cuid = %q, want proj-1", project.CUID)
	}
	if got := cloud.CountFor(http.MethodPost, "/api/projects/alice"); got != 0 {
		t.Fatalf("create requests = %d, want 0", got)
	}
}

// TestGetOrCreateProjectCreates 覆盖项目不存在：查询 404 → 新接口创建 → 重新查询。
func TestGetOrCreateProjectCreates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var created bool
	cloud.Handle(http.MethodGet, "/api/project/alice/demo", func(w http.ResponseWriter, _ *http.Request) {
		if !created {
			clienttest.WriteJSON(w, http.StatusNotFound, map[string]string{"code": "not_found"})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"cuid": "proj-1", "name": "demo"})
	})
	cloud.Handle(http.MethodPost, "/api/projects/alice", func(w http.ResponseWriter, _ *http.Request) {
		created = true
		clienttest.WriteJSON(w, http.StatusCreated, map[string]any{"cuid": "proj-1"})
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	project, err := svc.GetOrCreateProject(context.Background(), "alice", "demo", true)
	if err != nil {
		t.Fatalf("get or create project: %v", err)
	}
	if project.CUID != "proj-1" {
		t.Fatalf("project cuid = %q, want proj-1", project.CUID)
	}
	createReqs := cloud.RequestsFor(http.MethodPost, "/api/projects/alice")
	if len(createReqs) != 1 {
		t.Fatalf("create requests = %d, want 1", len(createReqs))
	}
	body := decodeJSONBody(t, createReqs[0].Body)
	if body["visibility"] != "PUBLIC" {
		t.Fatalf("visibility = %v, want PUBLIC", body["visibility"])
	}
	if got := cloud.CountFor(http.MethodGet, "/api/project/alice/demo"); got != 2 {
		t.Fatalf("project queries = %d, want 2 (before + after create)", got)
	}
}

// TestGetOrCreateProjectFallsBackToLegacy 覆盖旧后端：新接口 404 时调用旧接口，
// 409（已存在）视为成功。
func TestGetOrCreateProjectFallsBackToLegacy(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var created bool
	cloud.Handle(http.MethodGet, "/api/project/alice/demo", func(w http.ResponseWriter, _ *http.Request) {
		if !created {
			clienttest.WriteJSON(w, http.StatusNotFound, map[string]string{"code": "not_found"})
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"cuid": "proj-1", "name": "demo"})
	})
	cloud.ServeJSON(http.MethodPost, "/api/projects/alice", http.StatusNotFound, map[string]string{"code": "not_found"})
	cloud.Handle(http.MethodPost, "/api/project", func(w http.ResponseWriter, _ *http.Request) {
		created = true
		clienttest.WriteJSON(w, http.StatusConflict, map[string]string{"code": "conflict"})
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	if _, err := svc.GetOrCreateProject(context.Background(), "alice", "demo", false); err != nil {
		t.Fatalf("get or create project: %v", err)
	}
	legacyReqs := cloud.RequestsFor(http.MethodPost, "/api/project")
	if len(legacyReqs) != 1 {
		t.Fatalf("legacy create requests = %d, want 1", len(legacyReqs))
	}
	if body := decodeJSONBody(t, legacyReqs[0].Body); body["username"] != "alice" {
		t.Fatalf("legacy body username = %v, want alice", body["username"])
	}
}

// TestExperimentSummary 覆盖 Resume 摘要查询：createdAt 转成秒级时间戳作为查询下界。
func TestExperimentSummary(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-1", http.StatusOK, map[string]any{
		"log":    []map[string]any{{"key": "console", "step": 7}},
		"scalar": []map[string]any{{"key": "loss", "step": 12}},
		"media":  nil,
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	summary, err := svc.ExperimentSummary(context.Background(), "proj-1", "exp-1", "2026-01-02T03:04:05.000000Z")
	if err != nil {
		t.Fatalf("experiment summary: %v", err)
	}
	if len(summary.Scalar) != 1 || summary.Scalar[0].Key != "loss" || summary.Scalar[0].Step != 12 {
		t.Fatalf("scalar summary = %+v, want loss@12", summary.Scalar)
	}
	if len(summary.Log) != 1 || summary.Log[0].Step != 7 {
		t.Fatalf("log summary = %+v, want step 7", summary.Log)
	}

	reqs := cloud.RequestsFor(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-1")
	if len(reqs) != 1 {
		t.Fatalf("summary requests = %d, want 1", len(reqs))
	}
	wantCreatedAt := "1767323045" // 2026-01-02T03:04:05Z
	if got := reqs[0].Query["createdAt"]; len(got) != 1 || got[0] != wantCreatedAt {
		t.Fatalf("createdAt query = %v, want %s", got, wantCreatedAt)
	}
	if got := reqs[0].Query["all"]; len(got) != 1 || got[0] != "True" {
		t.Fatalf("all query = %v, want True (Python-compatible bool encoding)", got)
	}
}

// TestExperimentSummaryRejectsInvalidCreatedAt 覆盖 createdAt 非法：直接报错，
// 不回退到全历史扫描。
func TestExperimentSummaryRejectsInvalidCreatedAt(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	svc, _ := newTestClient(t, cloud, client.Options{})
	if _, err := svc.ExperimentSummary(context.Background(), "proj-1", "exp-1", ""); err == nil {
		t.Fatal("expected error for empty createdAt")
	}
	if got := cloud.CountFor(http.MethodGet, "/api/house/metrics/summaries/proj-1/exp-1"); got != 0 {
		t.Fatalf("summary requests = %d, want 0", got)
	}
}

// TestStopExperiment 覆盖终态上报：状态映射与 finishedAt 取请求时刻。
func TestStopExperiment(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPut, "/api/project/alice/demo/runs/exp-1/state", http.StatusOK, map[string]any{"ok": true})

	svc, _ := newTestClient(t, cloud, client.Options{})
	if err := svc.StopExperiment(context.Background(), StopExperimentParams{
		Username: "alice", Project: "demo", ExperimentID: "exp-1", State: RunStateCrashed,
	}); err != nil {
		t.Fatalf("stop experiment: %v", err)
	}
	reqs := cloud.RequestsFor(http.MethodPut, "/api/project/alice/demo/runs/exp-1/state")
	if len(reqs) != 1 {
		t.Fatalf("stop requests = %d, want 1", len(reqs))
	}
	body := decodeJSONBody(t, reqs[0].Body)
	if body["state"] != "CRASHED" {
		t.Fatalf("state = %v, want CRASHED", body["state"])
	}
	if body["from"] != "sdk" {
		t.Fatalf("from = %v, want sdk", body["from"])
	}
	if body["finishedAt"] != "2026-03-10T19:02:03.123400Z" {
		t.Fatalf("finishedAt = %v, want request-time UTC", body["finishedAt"])
	}
}

// TestUploadColumnsAndProfile 覆盖列定义与 Profile 上报的路径与请求体。
func TestUploadColumnsAndProfile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, "/api/projects/alice/demo/series", http.StatusOK, map[string]any{"ok": true})
	cloud.ServeJSON(http.MethodPut, "/api/project/alice/demo/runs/exp-1/profile", http.StatusOK, map[string]any{"ok": true})

	svc, _ := newTestClient(t, cloud, client.Options{})
	columns := []ColumnUpload{{Key: "loss", Type: "scalar"}, {Key: "img", Type: "image"}}
	if err := svc.UploadColumns(context.Background(), "alice", "demo", columns); err != nil {
		t.Fatalf("upload columns: %v", err)
	}
	seriesBody := decodeJSONBody(t, cloud.RequestsFor(http.MethodPost, "/api/projects/alice/demo/series")[0].Body)
	series, _ := seriesBody["series"].([]any)
	if len(series) != 2 {
		t.Fatalf("series = %v, want 2 columns", seriesBody["series"])
	}

	if err := svc.UploadProfile(context.Background(), UploadProfileParams{
		Username: "alice", Project: "demo", ExperimentID: "exp-1",
		Kind: ProfileConda, Content: "name: base",
	}); err != nil {
		t.Fatalf("upload profile: %v", err)
	}
	profileBody := decodeJSONBody(t, cloud.RequestsFor(http.MethodPut, "/api/project/alice/demo/runs/exp-1/profile")[0].Body)
	if profileBody["conda"] != "name: base" {
		t.Fatalf("profile body = %v, want conda content", profileBody)
	}
}

// TestRequestPresignedURLs 覆盖预签名申请：数量不匹配、空 URL 都是不可重试协议错误。
func TestRequestPresignedURLs(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	svc, _ := newTestClient(t, cloud, client.Options{})
	path := "/api/resources/presigned/put"

	cloud.ServeJSON(http.MethodPost, path, http.StatusOK, map[string]any{"urls": []string{"https://s3/a", "https://s3/b"}})
	urls, err := svc.RequestPresignedURLs(context.Background(), "exp-1", []string{"media/a.png", "media/b.png"})
	if err != nil {
		t.Fatalf("request presigned urls: %v", err)
	}
	if len(urls) != 2 || urls[1] != "https://s3/b" {
		t.Fatalf("urls = %v", urls)
	}

	cloud.ServeJSON(http.MethodPost, path, http.StatusOK, map[string]any{"urls": []string{"https://s3/a"}})
	_, err = svc.RequestPresignedURLs(context.Background(), "exp-1", []string{"a", "b"})
	if !client.IsInvalidResponse(err) || client.IsRetryable(err) {
		t.Fatalf("error = %v, want non-retryable invalid response", err)
	}

	cloud.ServeJSON(http.MethodPost, path, http.StatusOK, map[string]any{"urls": []string{"https://s3/a", ""}})
	_, err = svc.RequestPresignedURLs(context.Background(), "exp-1", []string{"a", "b"})
	if !client.IsInvalidResponse(err) {
		t.Fatalf("error = %v, want invalid response for empty url", err)
	}
}

// TestSaveFileEndpoints 覆盖 Save 小文件与前缀分片接口的请求体与响应校验。
func TestSaveFileEndpoints(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	svc, _ := newTestClient(t, cloud, client.Options{})

	cloud.ServeJSON(http.MethodPost, "/api/experiment/exp-1/files/prepare", http.StatusOK, map[string]any{
		"urls": []string{"https://s3/1"},
	})
	urls, err := svc.PrepareSaveFiles(context.Background(), "exp-1", []PrepareSaveFile{
		{Path: "files/a.txt", Size: 12, MD5: "abc", MIMEType: "text/plain"},
	})
	if err != nil {
		t.Fatalf("prepare save files: %v", err)
	}
	if len(urls) != 1 || urls[0] != "https://s3/1" {
		t.Fatalf("urls = %v", urls)
	}

	cloud.ServeJSON(http.MethodPost, "/api/experiment/exp-1/files/prepare", http.StatusOK, map[string]any{"urls": []string{}})
	if _, prepareErr := svc.PrepareSaveFiles(context.Background(), "exp-1", []PrepareSaveFile{{Path: "files/a.txt"}}); !client.IsInvalidResponse(prepareErr) {
		t.Fatalf("error = %v, want count mismatch as invalid response", prepareErr)
	}

	cloud.ServeJSON(http.MethodPost, "/api/experiment/exp-1/files/complete", http.StatusCreated, map[string]any{"ok": true})
	completedFiles := []CompleteSaveFile{{Path: "files/a.txt", State: MultipartUploaded}}
	if completeErr := svc.CompleteSaveFiles(context.Background(), "exp-1", completedFiles); completeErr != nil {
		t.Fatalf("complete save files: %v", completeErr)
	}

	cloud.ServeJSON(http.MethodPost, "/api/experiment/exp-1/files/prepare-multipart", http.StatusOK, map[string]any{
		"files": []map[string]any{{
			"uploadId": "up-1",
			"parts":    []map[string]any{{"partNumber": 1, "url": "https://s3/p1"}},
		}},
	})
	files, err := svc.PrepareMultipartSave(context.Background(), "exp-1", []PrepareMultipartSaveFile{
		{Path: "files/big.bin", Size: 100, MD5: "def", MIMEType: "application/octet-stream", Count: 1},
	})
	if err != nil {
		t.Fatalf("prepare multipart save: %v", err)
	}
	if len(files) != 1 || files[0].UploadID != "up-1" || files[0].Parts[0].PartNumber != 1 {
		t.Fatalf("multipart files = %+v", files)
	}
}

// TestCompleteMultipartParsesFileStates 覆盖 complete-multipart 的逐文件状态解析：
// 顶层数组、{"files": [...]} 与空响应体三种形态。
func TestCompleteMultipartParsesFileStates(t *testing.T) {
	path := "/api/experiment/exp-1/files/complete-multipart"
	completed := []CompletedMultipartSaveFile{{
		Path: "files/big.bin",
		Parts: []CompletedPart{
			{PartNumber: 1, ETag: "etag-1"},
		},
		UploadID: "up-1",
	}}

	t.Run("top level array", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		svc, _ := newTestClient(t, cloud, client.Options{})
		cloud.ServeJSON(http.MethodPost, path, http.StatusCreated, []map[string]string{
			{"path": "files/big.bin", "state": MultipartFailed},
		})
		states, err := svc.CompleteMultipartSave(context.Background(), "exp-1", completed)
		if err != nil {
			t.Fatalf("complete multipart: %v", err)
		}
		if len(states) != 1 || states[0].State != MultipartFailed {
			t.Fatalf("states = %+v, want one FAILED entry", states)
		}
	})

	t.Run("wrapped in files", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		svc, _ := newTestClient(t, cloud, client.Options{})
		cloud.ServeJSON(http.MethodPost, path, http.StatusCreated, map[string]any{
			"files": []map[string]string{{"path": "files/big.bin", "state": MultipartUploaded}},
		})
		states, err := svc.CompleteMultipartSave(context.Background(), "exp-1", completed)
		if err != nil {
			t.Fatalf("complete multipart: %v", err)
		}
		if len(states) != 1 || states[0].State != MultipartUploaded {
			t.Fatalf("states = %+v, want one UPLOADED entry", states)
		}
	})

	t.Run("empty body", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		svc, _ := newTestClient(t, cloud, client.Options{})
		cloud.Handle(http.MethodPost, path, func(w http.ResponseWriter, _ *http.Request) {
			w.WriteHeader(http.StatusCreated)
		})
		states, err := svc.CompleteMultipartSave(context.Background(), "exp-1", completed)
		if err != nil {
			t.Fatalf("complete multipart: %v", err)
		}
		if len(states) != 0 {
			t.Fatalf("states = %+v, want none for empty body", states)
		}
	})
}

// TestUploadMetricsAcceptedIsFailure 覆盖「假成功拦截」：House 返回 202 表示未写库，
// 必须判为不可重试失败。
func TestUploadMetricsAcceptedIsFailure(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.Handle(http.MethodPost, "/api/house/metrics", func(w http.ResponseWriter, _ *http.Request) {
		w.WriteHeader(http.StatusAccepted)
		_, _ = w.Write([]byte("session is invalid"))
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	err := svc.UploadMetrics(context.Background(), MetricScalar, "proj", "exp", []any{})
	if err == nil {
		t.Fatal("202 must not be treated as success")
	}
	if !client.IsHouseRejected(err) {
		t.Fatalf("error = %v, want house rejection", err)
	}
	if client.IsRetryable(err) {
		t.Fatal("house rejection must not be retryable")
	}
	if client.StatusCode(err) != http.StatusAccepted {
		t.Fatalf("status = %d, want 202", client.StatusCode(err))
	}
	if !strings.Contains(err.Error(), "session is invalid") {
		t.Fatalf("error = %v, want backend body preserved", err)
	}
}

// TestHeartbeatAcceptedIsFailure 覆盖心跳的同一规则（202 FlagId is invalid）。
func TestHeartbeatAcceptedIsFailure(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.Handle(http.MethodPost, "/api/house/experiments/exp-1/heartbeat", func(w http.ResponseWriter, _ *http.Request) {
		w.WriteHeader(http.StatusAccepted)
		_, _ = w.Write([]byte("FlagId is invalid"))
	})

	svc, _ := newTestClient(t, cloud, client.Options{})
	err := svc.SendHeartbeat(context.Background(), "exp-1")
	if !client.IsHouseRejected(err) || client.IsRetryable(err) {
		t.Fatalf("error = %v, want non-retryable house rejection", err)
	}
	if !strings.Contains(err.Error(), "FlagId is invalid") {
		t.Fatalf("error = %v, want backend body preserved", err)
	}
}

// TestFormatTimestamp 覆盖 createdAt/finishedAt 的 UTC 排版，与 Python
// Timestamp.ToDatetime().isoformat() + "Z" 对齐。
func TestFormatTimestamp(t *testing.T) {
	cases := []struct {
		name string
		in   time.Time
		want string
	}{
		{"with microseconds", clienttest.FixedNow, "2026-03-10T19:02:03.123400Z"},
		{"without microseconds", time.Date(2026, 3, 10, 19, 2, 3, 0, time.UTC), "2026-03-10T19:02:03Z"},
		{
			"converted to utc",
			time.Date(2026, 3, 11, 3, 2, 3, 0, time.FixedZone("CST", 8*3600)),
			"2026-03-10T19:02:03Z",
		},
	}
	for _, tc := range cases {
		if got := formatTimestamp(tc.in); got != tc.want {
			t.Fatalf("%s: formatTimestamp = %q, want %q", tc.name, got, tc.want)
		}
	}
}

// TestParseResponseTimestamp 覆盖摘要查询下界的时间戳解析。
func TestParseResponseTimestamp(t *testing.T) {
	want := clienttest.FixedNow.Unix()
	cases := []struct {
		in      string
		want    int64
		wantErr bool
	}{
		{in: "2026-03-10T19:02:03.123400Z", want: want},
		{in: "2026-03-10T19:02:03Z", want: want},
		{in: "2026-03-10T19:02:03", want: want},
		{in: "1773082923", want: 1773082923},
		{in: "1773082923123", want: 1773082923},
		{in: "", wantErr: true},
		{in: "not-a-timestamp", wantErr: true},
		{in: "123", wantErr: true},
	}
	for _, tc := range cases {
		got, err := parseResponseTimestamp(tc.in)
		if tc.wantErr {
			if err == nil {
				t.Fatalf("parseResponseTimestamp(%q) = %d, want error", tc.in, got)
			}
			continue
		}
		if err != nil {
			t.Fatalf("parseResponseTimestamp(%q): %v", tc.in, err)
		}
		if got != tc.want {
			t.Fatalf("parseResponseTimestamp(%q) = %d, want %d", tc.in, got, tc.want)
		}
	}
}
