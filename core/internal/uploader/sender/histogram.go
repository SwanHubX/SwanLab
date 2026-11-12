package sender

import (
	"context"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

const (
	// histogramChunk 是 /house/metrics 单次请求的直方图条数上限，与 Python 客户端一致。
	histogramChunk = 128
	// maxHistogramBuckets 是单条直方图的桶数上限。
	maxHistogramBuckets = 512

	// 本地校验失败时的确定性错误文案。
	histogramBlockedMessage      = "histogram record without value or timestamp"
	histogramInvalidShapeMessage = "histogram record has invalid bucket shape"
)

// HistogramSender 上报直方图指标：每批最多 128 条，data 只携带 values 与 bins。
type HistogramSender struct{ opts Options }

// NewHistogram 构造直方图 sender。
func NewHistogram(opts Options) *HistogramSender { return &HistogramSender{opts: opts} }

// Upload 实现 uploader.Sender。
func (s *HistogramSender) Upload(ctx context.Context, records []*recordv1.Record) uploader.Result {
	if len(records) == 0 {
		return uploader.Result{}
	}
	metrics := make([]map[string]any, 0, len(records))
	encodable := make([]*recordv1.Record, 0, len(records))
	blocked := make([]*recordv1.Record, 0)
	invalid := make([]*recordv1.Record, 0)
	for _, rec := range records {
		metric, status := encodeHistogram(rec.GetHistogram())
		switch status {
		case histogramOK:
			metrics = append(metrics, metric)
			encodable = append(encodable, rec)
		case histogramBadShape:
			invalid = append(invalid, rec)
		default:
			blocked = append(blocked, rec)
		}
	}

	failures := blockedFailures(histogramBlockedMessage, blocked)
	failures = append(failures, blockedFailures(histogramInvalidShapeMessage, invalid)...)
	terminal := append(append([]*recordv1.Record{}, blocked...), invalid...)
	if len(metrics) == 0 {
		return uploader.Result{Confirmed: records, Failures: failures}
	}

	var confirmed []*recordv1.Record
	for start := 0; start < len(encodable); start += histogramChunk {
		end := min(start+histogramChunk, len(encodable))
		batch := encodable[start:end]
		err := s.opts.API.UploadMetrics(ctx, api.MetricHistogram, s.opts.ProjectID, s.opts.ExperimentID,
			metrics[start:end])
		if err == nil {
			confirmed = append(confirmed, batch...)
			continue
		}
		failures = append(failures, uploader.Failure{Message: err.Error(), Records: batch, Err: err})
		if isRetryable(err) {
			return uploader.Result{
				Confirmed: append(confirmed, terminal...),
				Remaining: encodable[start:],
				Failures:  failures,
			}
		}
		// 确定性错误只终结当前分片：错误可能是记录级 400，后续分片仍须尝试，
		// 不能假定它们同样失败而把未发送的分片一并推进上传进度。
		confirmed = append(confirmed, batch...)
	}
	return uploader.Result{Confirmed: append(confirmed, terminal...), Failures: failures}
}

// histogramEncodeStatus 是直方图编码结果。
type histogramEncodeStatus int

const (
	histogramOK histogramEncodeStatus = iota
	histogramBlocked
	histogramBadShape
)

// encodeHistogram 编码单条直方图；缺 payload/时间戳与桶形状非法分别返回对应状态。
func encodeHistogram(record *datav1.HistogramRecord) (map[string]any, histogramEncodeStatus) {
	if record == nil || record.GetTimestamp() == nil || record.GetValue() == nil {
		return nil, histogramBlocked
	}
	createTime, ok := encodeTimestamp(record.GetTimestamp())
	if !ok {
		return nil, histogramBlocked
	}
	values := record.GetValue().GetValues()
	bins := record.GetValue().GetBins()
	if len(values) == 0 || len(values) > maxHistogramBuckets || len(bins) != len(values)+1 {
		return nil, histogramBadShape
	}
	return map[string]any{
		metricKeyField:   record.GetKey(),
		metricIndexField: record.GetStep(),
		metricDataField: map[string]any{
			"values": values,
			"bins":   bins,
		},
		createTimeKey: createTime,
	}, histogramOK
}
