package api

// jsonFilesKey 是 Save 协议请求体里的文件列表字段名，多处复用故提取为常量。
const jsonFilesKey = "files"

// ResumeMode 是实验恢复策略，取值与 proto ResumeMode 的语义一致。
type ResumeMode string

const (
	// ResumeNever 表示不复用已有实验。
	ResumeNever ResumeMode = "never"
	// ResumeAllow 表示实验存在则复用，不存在则新建。
	ResumeAllow ResumeMode = "allow"
	// ResumeMust 表示实验必须已存在，否则报错。
	ResumeMust ResumeMode = "must"
)

// RunState 是上报给云端的实验终态。
type RunState string

const (
	// RunStateFinished 表示正常结束。
	RunStateFinished RunState = "FINISHED"
	// RunStateCrashed 表示异常退出。
	RunStateCrashed RunState = "CRASHED"
	// RunStateAborted 表示用户主动中断。
	RunStateAborted RunState = "ABORTED"
)

// MetricKind 是 /house/metrics 的指标类型。
type MetricKind string

const (
	// MetricScalar 是标量指标。
	MetricScalar MetricKind = "scalar"
	// MetricLog 是控制台日志。
	MetricLog MetricKind = "log"
	// MetricMedia 是媒体指标。
	MetricMedia MetricKind = "media"
)

// ProfileKind 是实验 Profile 的四类元数据。
type ProfileKind string

const (
	// ProfileConda 是 conda 环境描述。
	ProfileConda ProfileKind = "conda"
	// ProfileRequirements 是 Python 依赖清单。
	ProfileRequirements ProfileKind = "requirements"
	// ProfileMetadata 是主机与环境元数据。
	ProfileMetadata ProfileKind = "metadata"
	// ProfileConfig 是实验配置。
	ProfileConfig ProfileKind = "config"
)

// 分片合并状态取值。
const (
	// MultipartUploaded 表示分片已合并成功。
	MultipartUploaded = "UPLOADED"
	// MultipartFailed 表示分片合并失败，需要重传。
	MultipartFailed = "FAILED"
)

// ---------------------------------------------------------------------------
// 响应体
// ---------------------------------------------------------------------------

// Experiment 是 create-or-resume 返回的实验信息。
type Experiment struct {
	CUID      string `json:"cuid"`
	Slug      string `json:"slug"`
	Name      string `json:"name"`
	CreatedAt string `json:"createdAt"`
}

// Project 是项目信息。Version 为 1 表示旧版（无多视图）项目。
type Project struct {
	CUID       string `json:"cuid"`
	Name       string `json:"name"`
	Path       string `json:"path"`
	Visibility string `json:"visibility"`
	Version    int    `json:"version"`
}

// ColumnStep 是 Resume 摘要中单列的步数水位。
type ColumnStep struct {
	Key  string `json:"key"`
	Step int64  `json:"step"`
}

// Summary 是 Resume 时拉取的实验摘要；三组列表分别对应日志、媒体与标量列。
type Summary struct {
	Log    []ColumnStep `json:"log"`
	Media  []ColumnStep `json:"media"`
	Scalar []ColumnStep `json:"scalar"`
}

// MultipartPart 是分片上传的单个分片信息。
type MultipartPart struct {
	PartNumber int    `json:"partNumber"`
	URL        string `json:"url"`
}

// MultipartFile 是 prepare-multipart 返回的单个文件分片信息。
type MultipartFile struct {
	UploadID string          `json:"uploadId"`
	Parts    []MultipartPart `json:"parts"`
}

// MultipartState 是 complete-multipart 响应中单个文件的合并结果。
type MultipartState struct {
	Path  string `json:"path"`
	State string `json:"state"`
}

// urlsResponse 是预签名与 Save prepare 共用的 {urls: [...]} 响应体。
type urlsResponse struct {
	URLs []string `json:"urls"`
}

// multipartFilesResponse 是 prepare-multipart 的 {files: [...]} 响应体。
type multipartFilesResponse struct {
	Files []MultipartFile `json:"files"`
}

// completedFilesResponse 是 complete-multipart 的 {"files": [...]} 响应体形态。
type completedFilesResponse struct {
	Files []MultipartState `json:"files"`
}

// ---------------------------------------------------------------------------
// 请求体
// ---------------------------------------------------------------------------

// CreateExperimentParams 是 create-or-resume 的入参。
//
// createdAt 不在这里暴露：后端按数据入库时间过滤查询，上报值必须是「发起请求的时刻」，
// 因此由 api 层在发送请求时取当前 UTC 时间（对齐 core_python/api/experiment.py）。
type CreateExperimentParams struct {
	Username    string
	Project     string
	Name        string
	Description string
	Color       string
	JobType     string
	Group       string
	Tags        []string
	RunID       string
	Resume      ResumeMode
}

// StopExperimentParams 是结束实验的入参。
type StopExperimentParams struct {
	Username     string
	Project      string
	ExperimentID string
	State        RunState
}

// UploadProfileParams 是 Profile 元数据上报的入参。Content 为文本（conda/requirements）
// 或已解析的结构（metadata/config）。
type UploadProfileParams struct {
	Username     string
	Project      string
	ExperimentID string
	Kind         ProfileKind
	Content      any
}

// ColumnUpload 是单列定义的上报形态（多视图之后的 series 协议）。
type ColumnUpload struct {
	Key         string `json:"key"`
	Type        string `json:"type"`
	SectionName string `json:"sectionName,omitempty"`
	XAxis       string `json:"xAxis,omitempty"`
	Hidden      bool   `json:"hidden,omitempty"`
}

// PrepareSaveFile 是单文件上传准备阶段（prepare）的单个文件描述。
type PrepareSaveFile struct {
	Path     string `json:"path"`
	Size     int64  `json:"size"`
	MD5      string `json:"md5"`
	MIMEType string `json:"mimeType"`
}

// CompleteSaveFile 是单文件上传完成阶段（complete）的单个文件状态。
type CompleteSaveFile struct {
	Path  string `json:"path"`
	State string `json:"state"`
}

// PrepareMultipartSaveFile 是分片上传准备阶段的单个文件描述。
type PrepareMultipartSaveFile struct {
	Path     string `json:"path"`
	Size     int64  `json:"size"`
	MD5      string `json:"md5"`
	MIMEType string `json:"mimeType"`
	Count    int    `json:"count"`
}

// CompletedPart 是 complete-multipart 上报的单个分片状态。
type CompletedPart struct {
	PartNumber int    `json:"partNumber"`
	ETag       string `json:"etag"`
}

// CompletedMultipartSaveFile 是 complete-multipart 上报的单个文件信息。
type CompletedMultipartSaveFile struct {
	Path     string          `json:"path"`
	Parts    []CompletedPart `json:"parts"`
	UploadID string          `json:"uploadId"`
}

// experimentBody 是 POST /project/{username}/{project}/experiment 的请求体。
//
// 空字段依靠 omitempty 剔除，等价于 Python 的 strip_none(strip_empty_str=True)。
type experimentBody struct {
	Name        string   `json:"name"`
	Description string   `json:"description,omitempty"`
	CreatedAt   string   `json:"createdAt"`
	Colors      []string `json:"colors"`
	Labels      []label  `json:"labels,omitempty"`
	Job         string   `json:"job,omitempty"`
	Cluster     string   `json:"cluster,omitempty"`
	CUID        string   `json:"cuid,omitempty"`
}

// label 是实验标签的请求体形态。
type label struct {
	Name string `json:"name"`
}
