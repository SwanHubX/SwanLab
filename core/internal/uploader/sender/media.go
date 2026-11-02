package sender

import (
	"bytes"
	"context"
	"io"
	"os"
	"path/filepath"
	"strconv"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/pkg/mime"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

const (
	// mediaMissingPayloadMessage 是记录缺少 media 载荷的文案。
	mediaMissingPayloadMessage = "media record without media payload"
	// mediaMissingTimestampMessage 是记录缺少时间戳的文案。
	mediaMissingTimestampMessage = "media record without timestamp"
	// mediaNoFilesMessage 是记录里没有任何可传文件的文案。
	mediaNoFilesMessage = "media record without files"
	// mediaNoPayloadMessage 是整批记录都没能收集到文件时的兜底文案。
	mediaNoPayloadMessage = "no media payload collected"
	// mediaSkipStoreReason 是 skip_store 下缺少 payload 的终结原因前缀。
	mediaSkipStoreReason = "payload missing with skip_store enabled: "
	// mediaFileNotFoundReason 是本地媒体文件不存在的终结原因前缀。
	mediaFileNotFoundReason = "media file not found: "
	// mediaUnsupportedTypeMessage 是未知媒体类型的终结文案前缀。
	mediaUnsupportedTypeMessage = "unsupported media type: "
)

// mediaMediums 是 MediaRecord.Type 到存储目录名（第二段路径）的映射。
var mediaMediums = map[columnv1.ColumnType]string{
	columnv1.ColumnType_COLUMN_TYPE_TEXT:     "text",
	columnv1.ColumnType_COLUMN_TYPE_IMAGE:    "image",
	columnv1.ColumnType_COLUMN_TYPE_AUDIO:    "audio",
	columnv1.ColumnType_COLUMN_TYPE_VIDEO:    "video",
	columnv1.ColumnType_COLUMN_TYPE_ECHARTS:  "echarts",
	columnv1.ColumnType_COLUMN_TYPE_OBJECT3D: "object3d",
	columnv1.ColumnType_COLUMN_TYPE_MOLECULE: "molecule",
	columnv1.ColumnType_COLUMN_TYPE_HTML:     "html",
}

// MediaSender 先把媒体文件直传对象存储，再提交引用这些文件的媒体指标。
type MediaSender struct{ opts Options }

// NewMedia 构造媒体 sender；缺少 API 或对象存储直传依赖时返回错误。
func NewMedia(opts Options) (*MediaSender, error) {
	if err := validateFileOptions(opts); err != nil {
		return nil, err
	}
	return &MediaSender{opts: opts}, nil
}

// mediaPayload 是单个媒体文件的上传载荷：内存字节或本地文件二选一。
type mediaPayload struct {
	remotePath  string
	trackerKey  string
	displayPath string
	size        int64
	contentType string
	taskIndex   int
	data        []byte // 非 nil 表示 payload 来自记录本身
	path        string // data 为 nil 时从本地文件流式读取
}

// mediaTask 是一条 media 记录的上传任务：metric 与 failure 互斥。
type mediaTask struct {
	index    int
	record   *recordv1.Record
	payloads []*mediaPayload
	metric   map[string]any
	failure  *uploader.Failure
}

// Upload 实现 uploader.Sender。
func (s *MediaSender) Upload(ctx context.Context, records []*recordv1.Record) uploader.Result {
	if len(records) == 0 {
		return uploader.Result{}
	}
	tasks, payloads, blocked, blockedFailures := s.collect(records)
	if len(payloads) == 0 {
		failures := blockedFailures
		if len(failures) == 0 {
			failures = []uploader.Failure{{Message: mediaNoPayloadMessage, Records: records}}
		}
		return uploader.Result{Confirmed: records, Failures: failures}
	}

	uploadable := make([]*recordv1.Record, 0, len(tasks)-len(blocked))
	for _, task := range tasks {
		if task.metric != nil {
			uploadable = append(uploadable, task.record)
		}
	}

	putErrors, err := s.putResources(ctx, payloads)
	if err != nil {
		return s.presignFailure(err, blocked, blockedFailures, uploadable, payloads)
	}

	confirmed := append([]*recordv1.Record{}, blocked...)
	failures := append([]uploader.Failure{}, blockedFailures...)
	var remaining []*recordv1.Record
	errorsByTask := make(map[int][]error, len(tasks))
	for i, putErr := range putErrors {
		if putErr != nil {
			errorsByTask[payloads[i].taskIndex] = append(errorsByTask[payloads[i].taskIndex], putErr)
		}
	}

	ready := make([]*mediaTask, 0, len(uploadable))
	for _, task := range tasks {
		if task.metric == nil {
			continue
		}
		taskErrors := errorsByTask[task.index]
		if len(taskErrors) == 0 {
			ready = append(ready, task)
			continue
		}
		terminal := false
		for _, taskErr := range taskErrors {
			failures = append(failures, uploader.Failure{
				Message: taskErr.Error(),
				Records: []*recordv1.Record{task.record},
				Err:     taskErr,
			})
			if !putErrorRetryable(taskErr) {
				terminal = true
			}
		}
		if terminal {
			confirmed = append(confirmed, task.record)
			for _, payload := range task.payloads {
				failFile(s.opts.Files.Tracker, payload.trackerKey)
			}
		} else {
			remaining = append(remaining, task.record)
		}
	}
	if len(ready) == 0 {
		return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
	}

	metrics := make([]map[string]any, 0, len(ready))
	readyRecords := make([]*recordv1.Record, 0, len(ready))
	for _, task := range ready {
		metrics = append(metrics, task.metric)
		readyRecords = append(readyRecords, task.record)
	}
	if err := s.opts.API.UploadMetrics(ctx, api.MetricMedia, s.opts.ProjectID, s.opts.ExperimentID, metrics); err != nil {
		failures = append(failures, uploader.Failure{Message: err.Error(), Records: readyRecords, Err: err})
		if isRetryable(err) {
			remaining = append(remaining, readyRecords...)
		} else {
			confirmed = append(confirmed, readyRecords...)
			for _, task := range ready {
				for _, payload := range task.payloads {
					failFile(s.opts.Files.Tracker, payload.trackerKey)
				}
			}
		}
		return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
	}
	for _, task := range ready {
		for _, payload := range task.payloads {
			finishFile(s.opts.Files.Tracker, payload.trackerKey)
		}
	}
	return uploader.Result{Confirmed: append(confirmed, readyRecords...), Remaining: remaining, Failures: failures}
}

// collect 把整批记录拆成上传任务：可传文件进 payloads，无文件的任务携带终结原因。
func (s *MediaSender) collect(records []*recordv1.Record) ([]*mediaTask, []*mediaPayload, []*recordv1.Record, []uploader.Failure) {
	tasks := make([]*mediaTask, 0, len(records))
	payloads := make([]*mediaPayload, 0, len(records))
	var blocked []*recordv1.Record
	var blockedFailures []uploader.Failure
	for _, rec := range records {
		task := &mediaTask{index: len(tasks), record: rec}
		tasks = append(tasks, task)
		s.collectTask(task)
		if task.failure != nil {
			blocked = append(blocked, rec)
			blockedFailures = append(blockedFailures, *task.failure)
			continue
		}
		payloads = append(payloads, task.payloads...)
	}
	return tasks, payloads, blocked, blockedFailures
}

// collectTask 收集单条记录的媒体载荷；无文件时写入 failure。
func (s *MediaSender) collectTask(task *mediaTask) {
	media := task.record.GetMedia()
	if media == nil {
		task.failure = blockedFailure(mediaMissingPayloadMessage, task.record)
		return
	}
	if media.GetTimestamp() == nil {
		task.failure = blockedFailure(mediaMissingTimestampMessage, task.record)
		return
	}
	medium, ok := mediaMediums[media.GetType()]
	if !ok {
		task.failure = blockedFailure(mediaUnsupportedTypeMessage+media.GetType().String(), task.record)
		return
	}
	createTime, ok := encodeTimestamp(media.GetTimestamp())
	if !ok {
		task.failure = blockedFailure(mediaMissingTimestampMessage, task.record)
		return
	}

	items := media.GetValue().GetItems()
	paths := make([]string, 0, len(items))
	more := make([]map[string]string, 0, len(items))
	reasons := make([]string, 0, len(items))
	for _, item := range items {
		remotePath := "media/" + medium + "/" + item.GetFilename()
		if item.Payload != nil {
			payload := &mediaPayload{
				remotePath:  remotePath,
				trackerKey:  remotePath + ":" + strconv.Itoa(len(item.Payload)),
				displayPath: remotePath,
				size:        int64(len(item.Payload)),
				contentType: mime.Type(item.GetFilename()),
				taskIndex:   task.index,
				data:        item.Payload,
			}
			task.payloads = append(task.payloads, payload)
			paths = append(paths, remotePath)
			if item.GetCaption() != "" {
				more = append(more, map[string]string{"caption": item.GetCaption()})
			}
			continue
		}
		if s.opts.Files.SkipStore {
			reasons = append(reasons, mediaSkipStoreReason+remotePath)
			continue
		}
		localPath := filepath.Join(s.opts.Files.MediaDir, medium, item.GetFilename())
		info, statErr := os.Stat(localPath)
		if statErr != nil || info.IsDir() {
			reasons = append(reasons, mediaFileNotFoundReason+localPath)
			continue
		}
		payload := &mediaPayload{
			remotePath:  remotePath,
			trackerKey:  remotePath + ":" + strconv.FormatInt(info.Size(), 10),
			displayPath: filepath.ToSlash(localPath),
			size:        info.Size(),
			contentType: mime.Type(localPath),
			taskIndex:   task.index,
			path:        localPath,
		}
		task.payloads = append(task.payloads, payload)
		paths = append(paths, remotePath)
		if item.GetCaption() != "" {
			more = append(more, map[string]string{"caption": item.GetCaption()})
		}
	}

	if len(paths) == 0 {
		reason := mediaNoFilesMessage
		if len(reasons) > 0 {
			reason = reasons[0]
		}
		task.failure = blockedFailure(reason, task.record)
		return
	}
	task.metric = map[string]any{
		"key":         media.GetKey(),
		"index":       media.GetStep(),
		"data":        paths,
		"more":        more,
		createTimeKey: createTime,
	}
}

// putResources 申请预签名 URL 并逐文件直传；返回按 payload 下标对齐的逐文件错误。
// 第二个返回值是整批失败（presign 不可用），此时不再尝试任何 PUT。
func (s *MediaSender) putResources(ctx context.Context, payloads []*mediaPayload) ([]error, error) {
	paths := make([]string, len(payloads))
	for i, payload := range payloads {
		paths[i] = payload.remotePath
	}
	for _, payload := range payloads {
		trackFile(s.opts.Files.Tracker, payload.trackerKey, payload.displayPath, payload.size)
	}

	urls, err := s.opts.API.RequestPresignedURLs(ctx, s.opts.ExperimentID, paths)
	if err != nil {
		return nil, err
	}
	errs := runConcurrently(ctx, mediaSlots, len(payloads), func(i int) error {
		return s.putPayload(ctx, payloads[i], urls[i])
	})
	return errs, nil
}

// putPayload 上传单个媒体文件，并在读取过程中上报字节进度；失败时回滚进度。
func (s *MediaSender) putPayload(ctx context.Context, payload *mediaPayload, url string) error {
	body, err := payload.open()
	if err != nil {
		return err
	}
	defer func() { _ = body.Close() }()

	reader := io.Reader(body)
	if tracker := s.opts.Files.Tracker; tracker != nil && payload.size > 0 {
		var current int64
		reader = &progressReader{reader: reader, onRead: func(n int64) {
			current += n
			tracker.UpdateFileProgress(payload.trackerKey, payload.trackerKey, current)
		}}
	}
	if _, err := s.opts.Files.Putter.PutObject(ctx, url, reader, payload.size, payload.contentType); err != nil {
		resetFileProgress(s.opts.Files.Tracker, payload.trackerKey, payload.trackerKey)
		return err
	}
	return nil
}

// open 返回载荷对应的读取器：内存字节或本地文件。
func (p *mediaPayload) open() (io.ReadCloser, error) {
	if p.data != nil {
		return io.NopCloser(bytes.NewReader(p.data)), nil
	}
	return os.Open(p.path)
}

// presignFailure 处理整批 presign 失败：可重试保留记录，否则终结并关闭文件行。
func (s *MediaSender) presignFailure(
	err error,
	blocked []*recordv1.Record,
	blockedFailures []uploader.Failure,
	uploadable []*recordv1.Record,
	payloads []*mediaPayload,
) uploader.Result {
	failures := append([]uploader.Failure{}, blockedFailures...)
	failures = append(failures, uploader.Failure{Message: err.Error(), Records: uploadable, Err: err})
	if isRetryable(err) {
		return uploader.Result{Confirmed: blocked, Remaining: uploadable, Failures: failures}
	}
	for _, payload := range payloads {
		failFile(s.opts.Files.Tracker, payload.trackerKey)
	}
	return uploader.Result{Confirmed: append(blocked, uploadable...), Failures: failures}
}
