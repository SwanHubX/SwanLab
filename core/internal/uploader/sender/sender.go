// Package sender 把 uploader 分发的整批记录编码成云端请求，并给出三态结果。
//
// 三态判定在这里完成：sender 持有 api 与 client.IsRetryable，Track 只做守恒校验
// 与记账。Confirmed 表示已终结（含 4xx、House 202、编码失败等确定性失败），
// Remaining 表示可重试。
package sender

import (
	"context"
	"errors"
	"fmt"
	"math"

	"google.golang.org/protobuf/encoding/protojson"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

const (
	// columnChunk 是 /series 单次请求的列数上限，对齐 Python upload_column。
	columnChunk = 3000
	// columnTypeSystem 是系统列在云端的类型标签，优先于真实类型。
	columnTypeSystem = "SYSTEM"
	// scalarBlockedMessage 对应 Python 的 blocked 分支文案。
	scalarBlockedMessage = "scalar record without value or timestamp"
	// logBlockedMessage 对应 Python 的 blocked 分支文案。
	logBlockedMessage = "log record without line or timestamp"
	// logUnknownLevelMessage 是 Go 侧新增的确定性失败：Python 遇到未知级别会崩溃。
	logUnknownLevelMessage = "log record with unknown level"
)

// columnTypeNames 是 ColumnType 到云端类型标签的映射，对齐 adapter.column。
var columnTypeNames = map[columnv1.ColumnType]string{
	columnv1.ColumnType_COLUMN_TYPE_SCALAR:   "FLOAT",
	columnv1.ColumnType_COLUMN_TYPE_IMAGE:    "IMAGE",
	columnv1.ColumnType_COLUMN_TYPE_AUDIO:    "AUDIO",
	columnv1.ColumnType_COLUMN_TYPE_VIDEO:    "VIDEO",
	columnv1.ColumnType_COLUMN_TYPE_TEXT:     "TEXT",
	columnv1.ColumnType_COLUMN_TYPE_ECHARTS:  "ECHARTS",
	columnv1.ColumnType_COLUMN_TYPE_OBJECT3D: "OBJECT3D",
	columnv1.ColumnType_COLUMN_TYPE_MOLECULE: "MOLECULE",
	columnv1.ColumnType_COLUMN_TYPE_HTML:     "HTML",
}

// logLevelNames 是 LogLevel 到云端日志级别的映射，对齐 adapter.level。
var logLevelNames = map[terminalv1.LogLevel]string{
	terminalv1.LogLevel_LOG_LEVEL_INFO:  "INFO",
	terminalv1.LogLevel_LOG_LEVEL_ERROR: "ERROR",
}

// Options 是三条 sender 共用的身份信息。
type Options struct {
	API          *api.Client // 必填
	Username     string
	Project      string
	ProjectID    string
	ExperimentID string
}

// ColumnSender 上报列定义：先全量编码，再按 3000 条切分提交。
type ColumnSender struct{ opts Options }

// NewColumn 构造列 sender。
func NewColumn(opts Options) *ColumnSender { return &ColumnSender{opts: opts} }

// Upload 实现 uploader.Sender。
func (s *ColumnSender) Upload(ctx context.Context, records []*recordv1.Record) uploader.Result {
	if len(records) == 0 {
		return uploader.Result{}
	}
	columns := make([]api.ColumnUpload, 0, len(records))
	for _, rec := range records {
		column, err := encodeColumn(rec.GetColumn())
		if err != nil {
			return uploader.Result{
				Confirmed: records,
				Failures: []uploader.Failure{{
					Message: "failed to encode column: " + err.Error(),
					Records: records,
					Err:     err,
				}},
			}
		}
		columns = append(columns, column)
	}

	var confirmed []*recordv1.Record
	var failures []uploader.Failure
	for start := 0; start < len(records); start += columnChunk {
		end := min(start+columnChunk, len(records))
		batch := records[start:end]
		if err := s.opts.API.UploadColumns(ctx, s.opts.Username, s.opts.Project, columns[start:end]); err != nil {
			failures = append(failures, uploader.Failure{Message: err.Error(), Records: batch, Err: err})
			if client.IsRetryable(err) {
				return uploader.Result{
					Confirmed: confirmed,
					Remaining: records[start:],
					Failures:  failures,
				}
			}
			// 确定性错误：剩余分片同样会失败，一并终结不再重复提交。
			return uploader.Result{
				Confirmed: append(confirmed, records[start:]...),
				Failures:  failures,
			}
		}
		confirmed = append(confirmed, batch...)
	}
	return uploader.Result{Confirmed: confirmed, Failures: failures}
}

// encodeColumn 把列记录转换成云端 DTO；未知类型返回错误。
func encodeColumn(record *columnv1.ColumnRecord) (api.ColumnUpload, error) {
	if record == nil {
		return api.ColumnUpload{}, errors.New("missing column payload")
	}
	column := api.ColumnUpload{Key: record.GetColumnKey()}
	if record.GetColumnClass() == columnv1.ColumnClass_COLUMN_CLASS_SYSTEM {
		column.Type = columnTypeSystem
	} else {
		typ, ok := columnTypeNames[record.GetColumnType()]
		if !ok {
			return api.ColumnUpload{}, fmt.Errorf("unsupported column type %s", record.GetColumnType())
		}
		column.Type = typ
	}
	if xAxis := record.GetXAxis(); xAxis != "" {
		column.XAxis = xAxis
	}
	if section := record.GetSectionName(); section != "" {
		column.SectionName = section
	}
	if record.GetHidden() {
		column.Hidden = true
	}
	return column, nil
}

// ScalarSender 上报标量指标，整批一次请求。
type ScalarSender struct{ opts Options }

// NewScalar 构造标量 sender。
func NewScalar(opts Options) *ScalarSender { return &ScalarSender{opts: opts} }

// Upload 实现 uploader.Sender。
func (s *ScalarSender) Upload(ctx context.Context, records []*recordv1.Record) uploader.Result {
	if len(records) == 0 {
		return uploader.Result{}
	}
	metrics := make([]map[string]any, 0, len(records))
	encodable := make([]*recordv1.Record, 0, len(records))
	blocked := make([]*recordv1.Record, 0)
	for _, rec := range records {
		metric, ok := encodeScalar(rec.GetScalar())
		if !ok {
			blocked = append(blocked, rec)
			continue
		}
		metrics = append(metrics, metric)
		encodable = append(encodable, rec)
	}
	if len(metrics) == 0 {
		return uploader.Result{
			Confirmed: records,
			Failures:  blockedFailures(scalarBlockedMessage, blocked),
		}
	}

	err := s.opts.API.UploadMetrics(ctx, api.MetricScalar, s.opts.ProjectID, s.opts.ExperimentID, metrics)
	if err == nil {
		return uploader.Result{
			Confirmed: append(encodable, blocked...),
			Failures:  blockedFailures(scalarBlockedMessage, blocked),
		}
	}

	failures := make([]uploader.Failure, 0, 1+len(blocked))
	failures = append(failures, uploader.Failure{Message: err.Error(), Records: encodable, Err: err})
	failures = append(failures, blockedFailures(scalarBlockedMessage, blocked)...)
	if client.IsRetryable(err) {
		return uploader.Result{Confirmed: blocked, Remaining: encodable, Failures: failures}
	}
	return uploader.Result{Confirmed: append(encodable, blocked...), Failures: failures}
}

// encodeScalar 编码单条标量；缺 value 或时间戳时返回 false。
func encodeScalar(record *datav1.ScalarRecord) (map[string]any, bool) {
	if record == nil || record.GetTimestamp() == nil || record.GetValue() == nil {
		return nil, false
	}
	createTime, ok := encodeTimestamp(record.GetTimestamp())
	if !ok {
		return nil, false
	}
	value := record.GetValue().GetNumber()
	var data any = value
	if math.IsNaN(value) || math.IsInf(value, 0) {
		data = "nan"
	}
	return map[string]any{
		"key":         record.GetKey(),
		"index":       record.GetStep(),
		"data":        data,
		"create_time": createTime,
	}, true
}

// LogSender 上报控制台日志，整批一次请求。
type LogSender struct{ opts Options }

// NewLog 构造日志 sender。
func NewLog(opts Options) *LogSender { return &LogSender{opts: opts} }

// Upload 实现 uploader.Sender。
func (s *LogSender) Upload(ctx context.Context, records []*recordv1.Record) uploader.Result {
	if len(records) == 0 {
		return uploader.Result{}
	}
	metrics := make([]map[string]any, 0, len(records))
	encodable := make([]*recordv1.Record, 0, len(records))
	blocked := make([]*recordv1.Record, 0)
	unknownLevel := make([]*recordv1.Record, 0)
	for _, rec := range records {
		metric, status := encodeLog(rec.GetLog())
		switch status {
		case logOK:
			metrics = append(metrics, metric)
			encodable = append(encodable, rec)
		case logBadLevel:
			unknownLevel = append(unknownLevel, rec)
		default:
			blocked = append(blocked, rec)
		}
	}

	terminal := append(blocked, unknownLevel...)
	blockedOnly := blockedFailures(logBlockedMessage, blocked)
	blockedOnly = append(blockedOnly, blockedFailures(logUnknownLevelMessage, unknownLevel)...)
	if len(metrics) == 0 {
		return uploader.Result{Confirmed: records, Failures: blockedOnly}
	}

	err := s.opts.API.UploadMetrics(ctx, api.MetricLog, s.opts.ProjectID, s.opts.ExperimentID, metrics)
	if err == nil {
		return uploader.Result{
			Confirmed: append(encodable, terminal...),
			Failures:  blockedOnly,
		}
	}

	failures := make([]uploader.Failure, 0, 1+len(blockedOnly))
	failures = append(failures, uploader.Failure{Message: err.Error(), Records: encodable, Err: err})
	failures = append(failures, blockedOnly...)
	if client.IsRetryable(err) {
		return uploader.Result{Confirmed: terminal, Remaining: encodable, Failures: failures}
	}
	return uploader.Result{Confirmed: append(encodable, terminal...), Failures: failures}
}

// logEncodeStatus 是日志编码结果。
type logEncodeStatus int

const (
	logOK logEncodeStatus = iota
	logBlocked
	logBadLevel
)

// encodeLog 编码单条日志；缺 payload/时间戳或级别未知时返回对应状态。
func encodeLog(record *terminalv1.LogRecord) (map[string]any, logEncodeStatus) {
	if record == nil || record.GetTimestamp() == nil {
		return nil, logBlocked
	}
	level, ok := logLevelNames[record.GetLevel()]
	if !ok {
		return nil, logBadLevel
	}
	createTime, ok := encodeTimestamp(record.GetTimestamp())
	if !ok {
		return nil, logBlocked
	}
	return map[string]any{
		"level":       level,
		"epoch":       record.GetEpoch(),
		"message":     record.GetLine(),
		"create_time": createTime,
	}, logOK
}

// encodeTimestamp 按 Python Timestamp.ToJsonString 的排版编码时间戳：
// 复用 protojson 保证 0/3/6/9 位小数规则一致，再去掉外层引号。
func encodeTimestamp(ts *timestamppb.Timestamp) (string, bool) {
	raw, err := protojson.Marshal(ts)
	if err != nil || len(raw) < 2 || raw[0] != '"' {
		return "", false
	}
	return string(raw[1 : len(raw)-1]), true
}

// blockedFailures 把一组本地校验失败的记录包装成 Failure；空集合返回 nil。
func blockedFailures(message string, records []*recordv1.Record) []uploader.Failure {
	if len(records) == 0 {
		return nil
	}
	return []uploader.Failure{{Message: message, Records: records}}
}
