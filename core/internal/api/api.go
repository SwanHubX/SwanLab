// Package api 封装 SwanLab 云端的 REST 业务协议：端点路径、请求体组装与响应校验
// （House 只认 200、实验创建区分 200/201、预签名数量校验）。
//
// 传输（连接池、鉴权、重试、错误模型）在 internal/client；协议模型在 model.go。
// 本包不感知会话、记录与上传管线，只把一次业务调用表达成成功或可分类的错误。
package api

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"net/http"
	"net/url"
	"strconv"
	"strings"
	"time"

	"github.com/swanhubx/swanlab/core/internal/client"
)

// 秒级 Unix 时间戳的合法区间，对齐 Python utils/time.TIMESTAMP_S_MIN/MAX。
const (
	timestampSecondsMin int64 = 1_000_000_000
	timestampSecondsMax int64 = 9_999_999_999
)

// Client 是云端 REST 端点的调用入口，绑定一个已认证的传输客户端。
//
// 并发安全：Client 自身无可变状态，可被同账号的多个会话共享。
type Client struct {
	http *client.Client
}

// New 用已认证的传输客户端构造 api 层。
func New(httpClient *client.Client) *Client {
	return &Client{http: httpClient}
}

// Username 返回鉴权得到的用户名，便于调用方省略显式的 owner 参数。
func (c *Client) Username() string { return c.http.Username() }

// ---------------------------------------------------------------------------
// 项目
// ---------------------------------------------------------------------------

// GetOrCreateProject 获取项目，不存在时创建，返回最终的项目信息。
//
// 流程：先查询，已存在直接返回（避免只写权限用户触发创建接口的 403）；404 时调用
// 新接口创建，新接口 404（旧后端）回退旧接口且 409 视为已存在；最后重新查询。
func (c *Client) GetOrCreateProject(ctx context.Context, username, name string, public bool) (*Project, error) {
	owner := username
	if owner == "" {
		owner = c.Username()
	}
	queryPath := projectPath(owner, name)

	project, err := c.getProject(ctx, queryPath)
	if err == nil {
		return project, nil
	}
	if !client.IsNotFound(err) {
		return nil, err
	}

	visibility := "PRIVATE"
	if public {
		visibility = "PUBLIC"
	}
	_, createErr := c.http.Do(ctx, client.Request{
		Method: http.MethodPost,
		Path:   "/projects/" + pathSegment(owner),
		Body:   map[string]string{"name": name, "visibility": visibility},
	})
	if createErr != nil {
		if !client.IsNotFound(createErr) {
			return nil, createErr
		}
		// 旧接口：已存在返回 409，同样视为成功。
		_, legacyErr := c.http.Do(ctx, client.Request{
			Method: http.MethodPost,
			Path:   "/project",
			Body: map[string]string{
				"name":       name,
				"visibility": visibility,
				"username":   owner,
			},
		})
		if legacyErr != nil && client.StatusCode(legacyErr) != http.StatusConflict {
			return nil, legacyErr
		}
	}
	return c.getProject(ctx, queryPath)
}

// getProject 查询单个项目。
func (c *Client) getProject(ctx context.Context, path string) (*Project, error) {
	res, err := c.http.Do(ctx, client.Request{Method: http.MethodGet, Path: path})
	if err != nil {
		return nil, err
	}
	var project Project
	if err := res.Decode(&project); err != nil {
		return nil, err
	}
	if project.CUID == "" {
		return nil, client.NewProtocolError(http.MethodGet, res.Endpoint, "invalid_project", "project response missing cuid", nil)
	}
	return &project, nil
}

// ---------------------------------------------------------------------------
// 实验
// ---------------------------------------------------------------------------

// CreateOrResumeExperiment 创建或恢复实验，返回实验信息以及是否为新建实验。
//
// resume=must 时先校验实验存在（404 即报错，避免静默新建出第二个实验）；
// 新建（201）时以本地上报的 createdAt 为准，因为旧后端 POST 响应不一定回传；
// 已存在（200）而响应缺 createdAt 时，回退到查询接口取真实创建时间。
func (c *Client) CreateOrResumeExperiment(
	ctx context.Context,
	params CreateExperimentParams,
) (*Experiment, bool, error) {
	if params.Resume == ResumeMust {
		if params.RunID == "" {
			return nil, false, errors.New("api: run_id is required when resume is \"must\"")
		}
		runPath := projectPath(params.Username, params.Project) + "/runs/" + pathSegment(params.RunID)
		if _, err := c.http.Do(ctx, client.Request{Method: http.MethodGet, Path: runPath}); err != nil {
			if client.IsNotFound(err) {
				return nil, false, fmt.Errorf(
					"experiment %s does not exist in project %s: %w", params.RunID, params.Project, err)
			}
			return nil, false, err
		}
	}

	reportedAt := formatTimestamp(c.http.Now())
	columns := []string{params.Color, params.Color}
	labels := make([]label, 0, len(params.Tags))
	for _, tag := range params.Tags {
		labels = append(labels, label{Name: tag})
	}
	body := experimentBody{
		Name:        params.Name,
		Description: params.Description,
		CreatedAt:   reportedAt,
		Colors:      columns,
		Labels:      labels,
		Job:         params.JobType,
		Cluster:     params.Group,
		CUID:        params.RunID,
	}

	res, err := c.http.Do(ctx, client.Request{
		Method: http.MethodPost,
		Path:   projectPath(params.Username, params.Project) + "/experiment",
		Body:   body,
	})
	if err != nil {
		return nil, false, err
	}
	if res.Status != http.StatusOK && res.Status != http.StatusCreated {
		return nil, false, client.NewProtocolError(
			res.Method, res.Endpoint, "unexpected_status",
			fmt.Sprintf("experiment create returned unexpected status %d", res.Status), nil)
	}

	var experiment Experiment
	if err := res.Decode(&experiment); err != nil {
		return nil, false, err
	}
	if experiment.CUID == "" {
		return nil, false, client.NewProtocolError(
			http.MethodPost, res.Endpoint, "invalid_experiment", "experiment response missing cuid", nil)
	}

	isNew := res.Status == http.StatusCreated
	if isNew {
		// 后端 POST 新建实验的响应不携带 createdAt，用上报时刻作为实验创建时间。
		experiment.CreatedAt = reportedAt
		return &experiment, true, nil
	}
	if experiment.CreatedAt == "" {
		// 旧后端 POST 响应未携带 createdAt，回退到查询接口。
		existing, getErr := c.getExperiment(ctx, params.Username, params.Project, experiment.CUID)
		if getErr != nil {
			return nil, false, getErr
		}
		experiment.CreatedAt = existing.CreatedAt
	}
	if experiment.CreatedAt == "" {
		return nil, false, client.NewProtocolError(
			http.MethodGet, res.Endpoint, "invalid_experiment",
			"backend did not return createdAt for experiment; please upgrade swanlab-server", nil)
	}
	return &experiment, false, nil
}

// getExperiment 按 run_id 查询实验。
func (c *Client) getExperiment(ctx context.Context, username, project, runID string) (*Experiment, error) {
	path := projectPath(username, project) + "/runs/" + pathSegment(runID)
	res, err := c.http.Do(ctx, client.Request{Method: http.MethodGet, Path: path})
	if err != nil {
		return nil, err
	}
	var experiment Experiment
	if err := res.Decode(&experiment); err != nil {
		return nil, err
	}
	return &experiment, nil
}

// ExperimentSummary 拉取实验摘要，用于 Resume 恢复列定义与步数水位。
//
// createdAt 作为数据入库时间的查询下界；缺失会让后端回退到全历史扫描，
// 因此解析失败直接报错，不做静默兜底。
func (c *Client) ExperimentSummary(
	ctx context.Context,
	projectID, experimentID, createdAt string,
) (*Summary, error) {
	seconds, err := parseResponseTimestamp(createdAt)
	if err != nil {
		return nil, fmt.Errorf("api: invalid experiment createdAt %q: %w", createdAt, err)
	}
	res, err := c.http.Do(ctx, client.Request{
		Method: http.MethodGet,
		Path:   fmt.Sprintf("/house/metrics/summaries/%s/%s", pathSegment(projectID), pathSegment(experimentID)),
		// "all" 取值与 Python requests 对 bool 的序列化保持一致（"True"）。
		Query: map[string]string{"all": "True", "createdAt": strconv.FormatInt(seconds, 10)},
	})
	if err != nil {
		return nil, err
	}
	var summary Summary
	if err := res.Decode(&summary); err != nil {
		return nil, err
	}
	return &summary, nil
}

// StopExperiment 上报实验终态。
//
// finishedAt 同样取发起请求时的当前 UTC 时间（后端按入库时间过滤查询，忽略历史时间戳）。
func (c *Client) StopExperiment(ctx context.Context, params StopExperimentParams) error {
	state := RunStateFinished
	switch params.State {
	case RunStateCrashed:
		state = RunStateCrashed
	case RunStateAborted:
		state = RunStateAborted
	}
	path := projectPath(params.Username, params.Project) + "/runs/" + pathSegment(params.ExperimentID) + "/state"
	_, err := c.http.Do(ctx, client.Request{
		Method: http.MethodPut,
		Path:   path,
		Body: map[string]string{
			"state":      string(state),
			"finishedAt": formatTimestamp(c.http.Now()),
			"from":       "sdk",
		},
	})
	return err
}

// SendHeartbeat 发送会话心跳，保持实验处于活跃状态。
//
// House 只把 200 视为成功：202 表示后端未写入（会话/flagId 校验失败），
// 当作成功会掩盖会话失效。Python 按 response.ok 判定，此处有意收紧。
func (c *Client) SendHeartbeat(ctx context.Context, experimentID string) error {
	res, err := c.http.Do(ctx, client.Request{
		Method: http.MethodPost,
		Path:   "/house/experiments/" + pathSegment(experimentID) + "/heartbeat",
		// flagId 是协议保留字段，Python 侧固定上报 "123456"。
		Body: map[string]string{"flagId": "123456"},
	})
	if err != nil {
		return err
	}
	if res.Status != http.StatusOK {
		return res.HouseRejected("heartbeat_rejected")
	}
	return nil
}

// ---------------------------------------------------------------------------
// 指标、列定义与 Profile
// ---------------------------------------------------------------------------

// UploadColumns 上报列定义（POST /projects/{username}/{project}/series）。
func (c *Client) UploadColumns(
	ctx context.Context,
	username, project string,
	columns []ColumnUpload,
) error {
	_, err := c.http.Do(ctx, client.Request{
		Method:  http.MethodPost,
		Path:    "/projects/" + pathSegment(username) + "/" + pathSegment(project) + "/series",
		Body:    map[string]any{"series": columns},
		Retries: client.NoRetries(),
	})
	return err
}

// UploadMetrics 把一批指标写入 House（POST /house/metrics）。
//
// metrics 由上层 sender 编码（标量/日志/媒体字段不同），本层只负责信封与结果判定：
// 仅 200 视为成功。空批次由 sender 过滤——metrics 是 any，判空需要反射。
func (c *Client) UploadMetrics(
	ctx context.Context,
	kind MetricKind,
	projectID, experimentID string,
	metrics any,
) error {
	res, err := c.http.Do(ctx, client.Request{
		Method: http.MethodPost,
		Path:   "/house/metrics",
		Body: map[string]any{
			"type":         string(kind),
			"projectId":    projectID,
			"experimentId": experimentID,
			"metrics":      metrics,
		},
		Retries: client.NoRetries(),
	})
	if err != nil {
		return err
	}
	if res.Status != http.StatusOK {
		return res.HouseRejected("session_rejected")
	}
	return nil
}

// UploadProfile 上报实验 Profile 元数据（四类合并到同一 PUT 接口的同名字段）。
func (c *Client) UploadProfile(ctx context.Context, params UploadProfileParams) error {
	path := projectPath(params.Username, params.Project) + "/runs/" + pathSegment(params.ExperimentID) + "/profile"
	_, err := c.http.Do(ctx, client.Request{
		Method:  http.MethodPut,
		Path:    path,
		Body:    map[string]any{string(params.Kind): params.Content},
		Retries: client.NoRetries(),
	})
	return err
}

// ---------------------------------------------------------------------------
// 文件上传（对象存储预签名与 Save 协议）
// ---------------------------------------------------------------------------

// RequestPresignedURLs 为一批对象存储路径申请预签名 PUT URL。
//
// 响应必须给出与 paths 等长、逐项非空的 urls 列表；数量或类型不符属于协议错误，
// 重试不会恢复，因此以不可重试错误终结该批记录。
func (c *Client) RequestPresignedURLs(ctx context.Context, experimentID string, paths []string) ([]string, error) {
	if len(paths) == 0 {
		return nil, nil
	}
	res, err := c.http.Do(ctx, client.Request{
		Method:  http.MethodPost,
		Path:    "/resources/presigned/put",
		Body:    map[string]any{"experimentId": experimentID, "paths": paths},
		Retries: client.NoRetries(),
	})
	if err != nil {
		return nil, err
	}
	var payload urlsResponse
	if err := res.Decode(&payload); err != nil {
		return nil, err
	}
	if err := validateURLs(res, "invalid_presigned_response", "presigned", payload.URLs, len(paths)); err != nil {
		return nil, err
	}
	return payload.URLs, nil
}

// PrepareSaveFiles 为小文件直传申请预签名 URL，返回顺序与入参一致。
func (c *Client) PrepareSaveFiles(
	ctx context.Context,
	experimentID string,
	files []PrepareSaveFile,
) ([]string, error) {
	res, err := c.http.Do(ctx, client.Request{
		Method:  http.MethodPost,
		Path:    "/experiment/" + pathSegment(experimentID) + "/files/prepare",
		Body:    map[string]any{jsonFilesKey: files},
		Retries: client.NoRetries(),
	})
	if err != nil {
		return nil, err
	}
	var payload urlsResponse
	if err := res.Decode(&payload); err != nil {
		return nil, err
	}
	if err := validateURLs(res, "invalid_prepare_response", "prepare", payload.URLs, len(files)); err != nil {
		return nil, err
	}
	return payload.URLs, nil
}

// CompleteSaveFiles 上报小文件直传结果。
func (c *Client) CompleteSaveFiles(
	ctx context.Context,
	experimentID string,
	files []CompleteSaveFile,
) error {
	_, err := c.http.Do(ctx, client.Request{
		Method:  http.MethodPost,
		Path:    "/experiment/" + pathSegment(experimentID) + "/files/complete",
		Body:    map[string]any{jsonFilesKey: files},
		Retries: client.NoRetries(),
	})
	return err
}

// PrepareMultipartSave 为大文件申请分片上传凭据，返回顺序与入参一致。
func (c *Client) PrepareMultipartSave(
	ctx context.Context,
	experimentID string,
	files []PrepareMultipartSaveFile,
) ([]MultipartFile, error) {
	res, err := c.http.Do(ctx, client.Request{
		Method:  http.MethodPost,
		Path:    "/experiment/" + pathSegment(experimentID) + "/files/prepare-multipart",
		Body:    map[string]any{jsonFilesKey: files},
		Retries: client.NoRetries(),
	})
	if err != nil {
		return nil, err
	}
	var payload multipartFilesResponse
	if err := res.Decode(&payload); err != nil {
		return nil, err
	}
	if err := validateMultipartFiles(res, payload.Files, files); err != nil {
		return nil, err
	}
	return payload.Files, nil
}

// CompleteMultipartSave 上报分片合并结果，并解析响应中的逐文件状态。
//
// 后端在合并失败时仍可能返回 201：若响应体给出逐文件状态列表则严格按状态分流；
// 响应体为空时返回 nil，由调用方按协议约定处理。
func (c *Client) CompleteMultipartSave(
	ctx context.Context,
	experimentID string,
	files []CompletedMultipartSaveFile,
) ([]MultipartState, error) {
	res, err := c.http.Do(ctx, client.Request{
		Method:  http.MethodPost,
		Path:    "/experiment/" + pathSegment(experimentID) + "/files/complete-multipart",
		Body:    map[string]any{jsonFilesKey: files},
		Retries: client.NoRetries(),
	})
	if err != nil {
		return nil, err
	}
	return decodeMultipartStates(res)
}

// decodeMultipartStates 兼容两种响应体形态：顶层数组，或 {"files": [...]}。
func decodeMultipartStates(res *client.Response) ([]MultipartState, error) {
	body := bytes.TrimSpace(res.Body)
	if len(body) == 0 {
		return nil, nil
	}
	if body[0] == '[' {
		var states []MultipartState
		if err := json.Unmarshal(body, &states); err != nil {
			return nil, client.NewProtocolError(res.Method, res.Endpoint, "invalid_complete_multipart_response", err.Error(), err)
		}
		return states, nil
	}
	var payload completedFilesResponse
	if err := json.Unmarshal(body, &payload); err != nil {
		return nil, client.NewProtocolError(res.Method, res.Endpoint, "invalid_complete_multipart_response", err.Error(), err)
	}
	return payload.Files, nil
}

// ---------------------------------------------------------------------------
// 内部工具
// ---------------------------------------------------------------------------

// validateURLs 校验预签名响应：与请求项等长，且逐项非空。
func validateURLs(res *client.Response, code, what string, urls []string, want int) error {
	if len(urls) != want {
		return client.NewProtocolError(res.Method, res.Endpoint, code,
			fmt.Sprintf("%s response returned %d urls, want %d", what, len(urls), want), nil)
	}
	for i, raw := range urls {
		if strings.TrimSpace(raw) == "" {
			return client.NewProtocolError(res.Method, res.Endpoint, code,
				fmt.Sprintf("%s response url[%d] is empty", what, i), nil)
		}
	}
	return nil
}

// validateMultipartFiles 校验 prepare-multipart 响应：文件数、uploadId 与分片编号。
func validateMultipartFiles(res *client.Response, files []MultipartFile, want []PrepareMultipartSaveFile) error {
	if len(files) != len(want) {
		return client.NewProtocolError(res.Method, res.Endpoint, "invalid_prepare_multipart_response",
			fmt.Sprintf("prepare-multipart response returned %d files, want %d", len(files), len(want)), nil)
	}
	for i, file := range files {
		if strings.TrimSpace(file.UploadID) == "" {
			return client.NewProtocolError(res.Method, res.Endpoint, "invalid_prepare_multipart_response",
				fmt.Sprintf("prepare-multipart file[%d] uploadId is empty", i), nil)
		}
		if err := validateParts(res, file.Parts, want[i].Count, i); err != nil {
			return err
		}
	}
	return nil
}

// validateParts 校验分片数与编号连续为 1..N，且每片 URL 非空。
func validateParts(res *client.Response, parts []MultipartPart, want, fileIndex int) error {
	invalid := func(message string) error {
		return client.NewProtocolError(res.Method, res.Endpoint, "invalid_prepare_multipart_response",
			fmt.Sprintf("prepare-multipart file[%d] %s", fileIndex, message), nil)
	}
	if len(parts) != want {
		return invalid(fmt.Sprintf("returned %d parts, want %d", len(parts), want))
	}
	seen := make([]bool, want+1)
	for _, part := range parts {
		if part.PartNumber < 1 || part.PartNumber > want || seen[part.PartNumber] {
			return invalid("part numbers are not consecutive from 1")
		}
		seen[part.PartNumber] = true
		if strings.TrimSpace(part.URL) == "" {
			return invalid(fmt.Sprintf("part[%d] url is empty", part.PartNumber))
		}
	}
	return nil
}

// projectPath 拼出 /project/{username}/{project} 路径。
func projectPath(username, project string) string {
	return "/project/" + pathSegment(username) + "/" + pathSegment(project)
}

// pathSegment 对单个路径段做转义，避免用户名/实验名中的特殊字符破坏 URL。
func pathSegment(value string) string { return url.PathEscape(value) }

// formatTimestamp 按 Python Timestamp→isoformat 的排版输出 UTC 时间：
// 微秒为 0 时不带小数部分，否则固定 6 位。
func formatTimestamp(at time.Time) string {
	at = at.UTC()
	base := at.Format("2006-01-02T15:04:05")
	if micro := at.Nanosecond() / 1000; micro > 0 {
		base += fmt.Sprintf(".%06d", micro)
	}
	return base + "Z"
}

// parseResponseTimestamp 把响应中的时间戳解析成秒级 Unix 时间戳。
//
// 支持数字（秒/毫秒）与 ISO 8601 字符串；无时区字符串按 UTC 解释。
// 对应 Python utils/time.parse_resp_timestamp_s 的语义与取值范围校验。
func parseResponseTimestamp(value string) (int64, error) {
	text := strings.TrimSpace(value)
	if text == "" {
		return 0, errors.New("timestamp value is required")
	}
	var seconds int64
	if digits := strings.TrimLeft(text, "-"); digits != "" && isDigits(digits) {
		parsed, err := strconv.ParseInt(text, 10, 64)
		if err != nil {
			return 0, fmt.Errorf("invalid numeric timestamp %q", text)
		}
		seconds = parsed
	} else {
		parsed, err := parseISO8601(text)
		if err != nil {
			return 0, err
		}
		seconds = parsed
	}
	if seconds <= 0 {
		return 0, fmt.Errorf("expected positive timestamp, got %d", seconds)
	}
	for seconds > timestampSecondsMax {
		seconds /= 1000
	}
	if seconds < timestampSecondsMin || seconds > timestampSecondsMax {
		return 0, fmt.Errorf("timestamp %d out of range [%d, %d]",
			seconds, timestampSecondsMin, timestampSecondsMax)
	}
	return seconds, nil
}

// parseISO8601 解析 ISO 8601 字符串并返回秒级 Unix 时间戳。
func parseISO8601(text string) (int64, error) {
	normalized := strings.Replace(text, "Z", "+00:00", 1)
	for _, layout := range []string{time.RFC3339Nano, "2006-01-02T15:04:05", "2006-01-02 15:04:05"} {
		if at, err := time.Parse(layout, normalized); err == nil {
			return at.Unix(), nil
		}
	}
	return 0, fmt.Errorf("invalid ISO 8601 timestamp %q", text)
}

// isDigits 判断字符串是否全为十进制数字。
func isDigits(text string) bool {
	if text == "" {
		return false
	}
	for _, r := range text {
		if r < '0' || r > '9' {
			return false
		}
	}
	return true
}
