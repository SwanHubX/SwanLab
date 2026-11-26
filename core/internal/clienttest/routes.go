package clienttest

import (
	"encoding/json"
	"net/http"
	"sync/atomic"
)

// ProjectSpec 描述假云端项目查询响应与元信息。
type ProjectSpec struct {
	Owner       string
	Name        string
	Cuid        string
	Path        string
	GroupOwner  string
	Version     int
	Experiments int
}

// ExperimentSpec 描述假云端实验创建/续传响应。
type ExperimentSpec struct {
	Owner     string
	Project   string
	Cuid      string
	Slug      string
	Status    int
	CreatedAt string
}

// RunSpec 描述假云端单次运行所需的默认路由集合。
type RunSpec struct {
	Owner       string
	Project     string
	ProjectCuid string
	ProjectPath string
	GroupOwner  string
	Version     int
	Experiments int
	ExpID       string
	Slug        string
	ExpStatus   int
	CreatedAt   string
	ServerUp    *atomic.Bool
	HouseUp     *atomic.Bool
	DownStatus  int
}

func orDefault(v, fallback string) string {
	if v != "" {
		return v
	}
	return fallback
}

func orDefaultInt(v, fallback int) int {
	if v != 0 {
		return v
	}
	return fallback
}

// ServeAuth 注册认证校验默认成功响应；username 为空时默认 "alice"。
func (c *Cloud) ServeAuth(username string) {
	username = orDefault(username, "alice")
	c.ServeJSON(http.MethodGet, "/api/auth/verify", http.StatusOK, map[string]any{
		"uid":       1,
		"username":  username,
		"createdAt": "2026-01-01T00:00:00Z",
	})
}

// ServeHealth 注册主站与指标服务的健康探针路由。
// 对应标志指针为 nil 时返回固定 200；非 nil 时按标志动态切换 200 与 downStatus（默认 503）。
func (c *Cloud) ServeHealth(serverUp, houseUp *atomic.Bool, downStatus int) {
	downStatus = orDefaultInt(downStatus, http.StatusServiceUnavailable)
	for _, ep := range []struct {
		path string
		up   *atomic.Bool
	}{
		{"/api/health", serverUp},
		{"/api/house/health", houseUp},
	} {
		flag := ep.up
		if flag == nil {
			c.ServeJSON(http.MethodGet, ep.path, http.StatusOK, map[string]any{"ok": true})
			continue
		}
		c.Handle(http.MethodGet, ep.path, func(w http.ResponseWriter, _ *http.Request) {
			if !flag.Load() {
				WriteJSON(w, downStatus, map[string]any{"ok": false})
				return
			}
			WriteJSON(w, http.StatusOK, map[string]any{"ok": true})
		})
	}
}

// ServeProject 注册项目查询响应。
func (c *Cloud) ServeProject(spec ProjectSpec) {
	owner := orDefault(spec.Owner, "alice")
	name := orDefault(spec.Name, "demo")
	groupOwner := orDefault(spec.GroupOwner, owner)
	c.ServeJSON(http.MethodGet, "/api/project/"+owner+"/"+name, http.StatusOK, map[string]any{
		"cuid":       orDefault(spec.Cuid, "proj-1"),
		"name":       name,
		"path":       orDefault(spec.Path, "/"+owner+"/"+name),
		"visibility": "PRIVATE",
		"version":    orDefaultInt(spec.Version, 2),
		"_count":     map[string]any{"experiments": orDefaultInt(spec.Experiments, 1)},
		"group":      map[string]any{"username": groupOwner},
	})
}

// ServeExperiment 注册实验创建/续传路由，回显请求中的 name。
func (c *Cloud) ServeExperiment(spec ExperimentSpec) {
	owner := orDefault(spec.Owner, "alice")
	project := orDefault(spec.Project, "demo")
	cuid := orDefault(spec.Cuid, "exp-1")
	slug := orDefault(spec.Slug, "run-1")
	status := orDefaultInt(spec.Status, http.StatusCreated)
	c.Handle(http.MethodPost, "/api/project/"+owner+"/"+project+"/experiment", func(w http.ResponseWriter, r *http.Request) {
		var body struct {
			Name string `json:"name"`
		}
		_ = json.NewDecoder(r.Body).Decode(&body)
		payload := map[string]any{
			"cuid": cuid,
			"slug": slug,
			"name": body.Name,
		}
		if spec.CreatedAt != "" {
			payload["createdAt"] = spec.CreatedAt
		}
		WriteJSON(w, status, payload)
	})
}

// ServeHeartbeat 注册实验心跳成功响应。
func (c *Cloud) ServeHeartbeat(experimentID string) {
	c.ServeJSON(http.MethodPost, "/api/house/experiments/"+orDefault(experimentID, "exp-1")+"/heartbeat",
		http.StatusOK, map[string]any{"ok": true})
}

// ServeSeries 注册列元数据上报成功响应。
func (c *Cloud) ServeSeries(owner, project string) {
	c.ServeJSON(http.MethodPost, "/api/projects/"+orDefault(owner, "alice")+"/"+orDefault(project, "demo")+"/series",
		http.StatusOK, map[string]any{"ok": true})
}

// ServeRunState 注册实验终态上报成功响应。
func (c *Cloud) ServeRunState(owner, project, experimentID string) {
	c.ServeJSON(http.MethodPut,
		"/api/project/"+orDefault(owner, "alice")+"/"+orDefault(project, "demo")+"/runs/"+orDefault(experimentID, "exp-1")+"/state",
		http.StatusOK, map[string]any{"ok": true})
}

// ServeSummary 注册实验指标摘要响应；summary 为 nil 时返回空标量摘要。
func (c *Cloud) ServeSummary(projectID, experimentID string, summary map[string]any) {
	if summary == nil {
		summary = map[string]any{"scalar": []any{}}
	}
	c.ServeJSON(http.MethodGet,
		"/api/house/metrics/summaries/"+orDefault(projectID, "proj-1")+"/"+orDefault(experimentID, "exp-1"),
		http.StatusOK, summary)
}

// ServeMetricsOK 注册 house 指标写入成功响应。
func (c *Cloud) ServeMetricsOK() {
	c.ServeJSON(http.MethodPost, "/api/house/metrics", http.StatusOK, map[string]any{"ok": true})
}

// ServeRun 一键注册单次运行所需的标准 SwanLab 后端路由集合；调用后仍可用 Handle/ServeJSON 覆写任意路由。
func (c *Cloud) ServeRun(spec RunSpec) {
	owner := orDefault(spec.Owner, "alice")
	project := orDefault(spec.Project, "demo")
	groupOwner := orDefault(spec.GroupOwner, owner)
	expID := orDefault(spec.ExpID, "exp-1")

	c.ServeAuth("")
	c.ServeHealth(spec.ServerUp, spec.HouseUp, spec.DownStatus)
	c.ServeProject(ProjectSpec{
		Owner:       owner,
		Name:        project,
		Cuid:        spec.ProjectCuid,
		Path:        spec.ProjectPath,
		GroupOwner:  groupOwner,
		Version:     spec.Version,
		Experiments: spec.Experiments,
	})
	c.ServeExperiment(ExperimentSpec{
		Owner:     groupOwner,
		Project:   project,
		Cuid:      expID,
		Slug:      spec.Slug,
		Status:    spec.ExpStatus,
		CreatedAt: spec.CreatedAt,
	})
	c.ServeHeartbeat(expID)
	c.ServeSeries(groupOwner, project)
	c.ServeRunState(groupOwner, project, expID)
	c.ServeMetricsOK()
}
