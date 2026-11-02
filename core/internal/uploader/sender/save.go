package sender

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"unicode/utf8"

	"gopkg.in/yaml.v3"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/fs"
	"github.com/swanhubx/swanlab/core/internal/pkg/mime"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

const (
	// saveNoPayloadMessage 是记录缺少 save 载荷的文案。
	saveNoPayloadMessage = "record without save payload"
	// saveCustomPayloadMessage 是 CUSTOM 记录携带 payload 的协议错误文案。
	saveCustomPayloadMessage = "CUSTOM save must not carry payload"
	// saveFileNotFoundPrefix 是 CUSTOM 源文件缺失的终结原因前缀。
	saveFileNotFoundPrefix = "save file not found: "
	// saveFileUnreadablePrefix 是源文件 stat 失败的终结原因前缀。
	saveFileUnreadablePrefix = "save file unreadable: "
	// saveFileTooLargePrefix 是超过单文件上限的终结原因前缀。
	saveFileTooLargePrefix = "save file exceeds size limit: "
	// saveMD5FailedPrefix 是 MD5 计算失败的终结原因前缀。
	saveMD5FailedPrefix = "save file md5 computation failed: "
	// saveUploadFailedPrefix 是单文件直传失败的诊断前缀。
	saveUploadFailedPrefix = "save file upload failed: "
	// saveMultipartFailedPrefix 是分片上传失败的诊断前缀。
	saveMultipartFailedPrefix = "multipart upload failed: "
	// saveMultipartStateFailedPrefix 是 complete-multipart 返回 FAILED 的诊断前缀。
	saveMultipartStateFailedPrefix = "multipart save failed: "
	// saveInternalUnavailablePrefix 是内部 save 内容不可用的终结原因前缀。
	saveInternalUnavailablePrefix = "internal save content unavailable: "
	// saveInternalUnparseablePrefix 是内部 save 解析失败的终结原因前缀。
	saveInternalUnparseablePrefix = "internal save payload unparseable: "

	// saveStateUploaded 是 complete 阶段上报的成功状态。
	saveStateUploaded = "UPLOADED"
	// saveStateFailed 是 complete 阶段上报的本地失败状态。
	saveStateFailed = "FAILED"
	// savePartContentType 是分片上传的内容类型。
	savePartContentType = "application/octet-stream"
)

// SaveSender 分发内部元文件与 CUSTOM 文件：前者进 Profile，后者直传对象存储。
type SaveSender struct{ opts Options }

// NewSave 构造 save sender；缺少 API/直传依赖或 Save 参数非正时返回错误。
func NewSave(opts Options) (*SaveSender, error) {
	if err := validateFileOptions(opts); err != nil {
		return nil, err
	}
	files := opts.Files
	switch {
	case files.SaveSize <= 0:
		return nil, fmt.Errorf("sender: save size must be positive, got %d", files.SaveSize)
	case files.SaveSplit <= 0:
		return nil, fmt.Errorf("sender: save split must be positive, got %d", files.SaveSplit)
	case files.SavePart <= 0:
		return nil, fmt.Errorf("sender: save part must be positive, got %d", files.SavePart)
	case files.SaveBatch <= 0:
		return nil, fmt.Errorf("sender: save batch must be positive, got %d", files.SaveBatch)
	}
	return &SaveSender{opts: opts}, nil
}

// pendingSave 是通过校验、等待上传的 CUSTOM 文件。
type pendingSave struct {
	record     *recordv1.Record
	sourcePath string
	size       int64
	name       string
}

// saveEntry 是 prepare 与 tracker 使用的文件描述。
type saveEntry struct {
	path       string
	sourcePath string
	size       int64
	md5        string
	mimeType   string
}

// saveTask 是单个 CUSTOM 文件的上传任务。
type saveTask struct {
	record     *recordv1.Record
	entry      saveEntry
	trackerKey string
}

// multipartReady 是分片已全部上传、等待合并的文件。
type multipartReady struct {
	task     *saveTask
	uploadID string
	parts    []api.CompletedPart
}

// Upload 实现 uploader.Sender。
func (s *SaveSender) Upload(ctx context.Context, records []*recordv1.Record) uploader.Result {
	if len(records) == 0 {
		return uploader.Result{}
	}
	confirmed := make([]*recordv1.Record, 0, len(records))
	var remaining []*recordv1.Record
	var failures []uploader.Failure
	custom := make([]*recordv1.Record, 0, len(records))
	for _, record := range records {
		save := record.GetSave()
		if save == nil {
			failures = append(failures, uploader.Failure{Message: saveNoPayloadMessage, Records: []*recordv1.Record{record}})
			confirmed = append(confirmed, record)
			continue
		}
		if save.GetType() == savev1.SaveType_SAVE_TYPE_CUSTOM {
			if save.Payload != nil {
				failures = append(failures, uploader.Failure{
					Message: saveCustomPayloadMessage,
					Records: []*recordv1.Record{record},
				})
				confirmed = append(confirmed, record)
				continue
			}
			custom = append(custom, record)
			continue
		}
		result := s.uploadInternal(ctx, record)
		confirmed = append(confirmed, result.Confirmed...)
		remaining = append(remaining, result.Remaining...)
		failures = append(failures, result.Failures...)
	}
	if len(custom) > 0 {
		result := s.uploadCustom(ctx, custom)
		confirmed = append(confirmed, result.Confirmed...)
		remaining = append(remaining, result.Remaining...)
		failures = append(failures, result.Failures...)
	}
	return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
}

// uploadInternal 读取并解码一条内部 save，合并到 Profile 上报。
func (s *SaveSender) uploadInternal(ctx context.Context, record *recordv1.Record) uploader.Result {
	save := record.GetSave()
	data, ok := s.readInternal(save)
	if !ok {
		return confirmedFailure(record, saveInternalUnavailablePrefix+save.GetSourcePath())
	}
	kind, content, ok := decodeInternal(save.GetType(), data)
	if !ok {
		return confirmedFailure(record, saveInternalUnparseablePrefix+save.GetType().String())
	}
	err := s.opts.API.UploadProfile(ctx, api.UploadProfileParams{
		Username:     s.opts.Username,
		Project:      s.opts.Project,
		ExperimentID: s.opts.ExperimentID,
		Kind:         kind,
		Content:      content,
	})
	if err != nil {
		failure := uploader.Failure{Message: err.Error(), Records: []*recordv1.Record{record}, Err: err}
		if isRetryable(err) {
			return uploader.Result{Remaining: []*recordv1.Record{record}, Failures: []uploader.Failure{failure}}
		}
		return uploader.Result{Confirmed: []*recordv1.Record{record}, Failures: []uploader.Failure{failure}}
	}
	return uploader.Result{Confirmed: []*recordv1.Record{record}}
}

// readInternal 返回内部 save 的字节内容：payload 优先，其次本地文件。
// skip_store 下没有本 run 的镜像，禁止回退到可能属于旧 run 的同名文件。
func (s *SaveSender) readInternal(save *savev1.SaveRecord) ([]byte, bool) {
	if save.Payload != nil {
		return save.Payload, true
	}
	if s.opts.Files.SkipStore {
		return nil, false
	}
	source := s.resolveSource(save)
	if source == "" {
		return nil, false
	}
	data, err := os.ReadFile(source)
	if err != nil {
		return nil, false
	}
	return data, true
}

// decodeInternal 按 SaveType 解码内部 save；返回 Profile 类型与上报内容。
func decodeInternal(saveType savev1.SaveType, data []byte) (api.ProfileKind, any, bool) {
	if !utf8.Valid(data) {
		return "", nil, false
	}
	switch saveType {
	case savev1.SaveType_SAVE_TYPE_METADATA:
		content, ok := decodeMetadata(data)
		if !ok {
			return "", nil, false
		}
		return api.ProfileMetadata, content, true
	case savev1.SaveType_SAVE_TYPE_CONFIG:
		var content map[string]any
		if err := yaml.Unmarshal(data, &content); err != nil || content == nil {
			return "", nil, false
		}
		return api.ProfileConfig, content, true
	case savev1.SaveType_SAVE_TYPE_REQUIREMENTS:
		if len(data) == 0 {
			return "", nil, false
		}
		return api.ProfileRequirements, string(data), true
	case savev1.SaveType_SAVE_TYPE_CONDA:
		if len(data) == 0 {
			return "", nil, false
		}
		return api.ProfileConda, string(data), true
	}
	return "", nil, false
}

// decodeMetadata 解码 metadata JSON：保留整数精度，且只接受单个 JSON 值。
func decodeMetadata(data []byte) (map[string]any, bool) {
	decoder := json.NewDecoder(bytes.NewReader(data))
	decoder.UseNumber()
	var content map[string]any
	if err := decoder.Decode(&content); err != nil || content == nil {
		return nil, false
	}
	if _, err := decoder.Token(); !errors.Is(err, io.EOF) {
		return nil, false
	}
	return content, true
}

// uploadCustom 解析、校验并分批上传 CUSTOM 文件。
func (s *SaveSender) uploadCustom(ctx context.Context, records []*recordv1.Record) uploader.Result {
	confirmed := make([]*recordv1.Record, 0, len(records))
	var failures []uploader.Failure
	pending := make([]pendingSave, 0, len(records))
	for _, record := range records {
		save := record.GetSave()
		source := s.resolveSource(save)
		if source == "" {
			failures = append(failures, uploader.Failure{
				Message: saveFileNotFoundPrefix + save.GetSourcePath(),
				Records: []*recordv1.Record{record},
			})
			confirmed = append(confirmed, record)
			continue
		}
		info, err := os.Stat(source)
		if err != nil {
			failures = append(failures, uploader.Failure{
				Message: saveFileUnreadablePrefix + source,
				Records: []*recordv1.Record{record},
				Err:     err,
			})
			confirmed = append(confirmed, record)
			continue
		}
		if info.Size() > s.opts.Files.SaveSize {
			failures = append(failures, uploader.Failure{
				Message: saveFileTooLargePrefix + source,
				Records: []*recordv1.Record{record},
			})
			confirmed = append(confirmed, record)
			continue
		}
		pending = append(pending, pendingSave{
			record:     record,
			sourcePath: source,
			size:       info.Size(),
			name:       save.GetName(),
		})
	}
	if len(pending) == 0 {
		return uploader.Result{Confirmed: confirmed, Failures: failures}
	}
	batch := s.opts.Files.SaveBatch
	var remaining []*recordv1.Record
	for start := 0; start < len(pending); start += batch {
		end := min(start+batch, len(pending))
		result := s.uploadBatch(ctx, pending[start:end])
		confirmed = append(confirmed, result.Confirmed...)
		remaining = append(remaining, result.Remaining...)
		failures = append(failures, result.Failures...)
	}
	return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
}

// uploadBatch 计算 MD5 后上传一批 CUSTOM 文件。
func (s *SaveSender) uploadBatch(ctx context.Context, pending []pendingSave) uploader.Result {
	paths := make([]string, len(pending))
	for i, item := range pending {
		paths[i] = item.sourcePath
	}
	digests, errs := computeMD5s(ctx, paths)
	confirmed := make([]*recordv1.Record, 0, len(pending))
	var remaining []*recordv1.Record
	var failures []uploader.Failure
	tasks := make([]*saveTask, 0, len(pending))
	for i, item := range pending {
		if errs[i] != nil {
			failure := uploader.Failure{
				Message: saveMD5FailedPrefix + item.sourcePath,
				Records: []*recordv1.Record{item.record},
				Err:     errs[i],
			}
			failures = append(failures, failure)
			if isRetryable(errs[i]) {
				remaining = append(remaining, item.record)
			} else {
				confirmed = append(confirmed, item.record)
			}
			continue
		}
		entry := saveEntry{
			path:       item.name,
			sourcePath: item.sourcePath,
			size:       item.size,
			md5:        digests[i],
			mimeType:   mime.Type(item.sourcePath),
		}
		tasks = append(tasks, &saveTask{record: item.record, entry: entry, trackerKey: saveTrackerKey(entry)})
	}
	if len(tasks) == 0 {
		return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
	}
	split := s.opts.Files.SaveSplit
	small := make([]*saveTask, 0, len(tasks))
	large := make([]*saveTask, 0)
	for _, task := range tasks {
		if task.entry.size < split {
			small = append(small, task)
		} else {
			large = append(large, task)
		}
	}
	result := uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
	if len(small) > 0 {
		result = mergeResult(result, s.uploadSmall(ctx, small))
	}
	if len(large) > 0 {
		result = mergeResult(result, s.uploadLarge(ctx, large))
	}
	return result
}

// uploadSmall 走 prepare → 并发 PUT → complete 的小文件路径。
func (s *SaveSender) uploadSmall(ctx context.Context, tasks []*saveTask) uploader.Result {
	items := make([]api.PrepareSaveFile, len(tasks))
	for i, task := range tasks {
		items[i] = api.PrepareSaveFile{
			Path:     task.entry.path,
			Size:     task.entry.size,
			MD5:      task.entry.md5,
			MIMEType: task.entry.mimeType,
		}
	}
	urls, err := s.opts.API.PrepareSaveFiles(ctx, s.opts.ExperimentID, items)
	if err != nil {
		return batchFailure(taskRecords(tasks), err)
	}
	for _, task := range tasks {
		trackFile(s.opts.Files.Tracker, task.trackerKey, task.entry.sourcePath, task.entry.size)
	}

	putErrors := runConcurrently(ctx, saveSlots, len(tasks), func(i int) error {
		return s.putSave(ctx, tasks[i], urls[i])
	})
	uploaded := make([]*saveTask, 0, len(tasks))
	storageFailed := make([]*saveTask, 0)
	localFailed := make([]*saveTask, 0)
	var failures []uploader.Failure
	for i, task := range tasks {
		if putErrors[i] == nil {
			uploaded = append(uploaded, task)
			continue
		}
		failures = append(failures, uploader.Failure{
			Message: saveUploadFailedPrefix + task.entry.path,
			Records: []*recordv1.Record{task.record},
			Err:     putErrors[i],
		})
		if putErrorRetryable(putErrors[i]) {
			storageFailed = append(storageFailed, task)
		} else {
			localFailed = append(localFailed, task)
		}
	}

	completeItems := make([]api.CompleteSaveFile, 0, len(uploaded)+len(localFailed))
	for _, task := range uploaded {
		completeItems = append(completeItems, api.CompleteSaveFile{Path: task.entry.path, State: saveStateUploaded})
	}
	for _, task := range localFailed {
		completeItems = append(completeItems, api.CompleteSaveFile{Path: task.entry.path, State: saveStateFailed})
	}
	if len(completeItems) > 0 {
		if err := s.opts.API.CompleteSaveFiles(ctx, s.opts.ExperimentID, completeItems); err != nil {
			return s.completeFailure(err, uploaded, storageFailed, localFailed, failures)
		}
	}
	for _, task := range uploaded {
		finishFile(s.opts.Files.Tracker, task.trackerKey)
	}
	for _, task := range localFailed {
		failFile(s.opts.Files.Tracker, task.trackerKey)
	}
	if len(uploaded) > 0 {
		console.Infof("Synced %d save file(s)", len(uploaded))
	}
	confirmed := append(taskRecords(uploaded), taskRecords(localFailed)...)
	return uploader.Result{Confirmed: confirmed, Remaining: taskRecords(storageFailed), Failures: failures}
}

// completeFailure 处理 complete 失败：可重试保留整批，否则终结并关闭进度行。
func (s *SaveSender) completeFailure(
	err error,
	uploaded, storageFailed, localFailed []*saveTask,
	failures []uploader.Failure,
) uploader.Result {
	pending := make([]*saveTask, 0, len(uploaded)+len(storageFailed))
	pending = append(pending, uploaded...)
	pending = append(pending, storageFailed...)
	pendingRecords := taskRecords(pending)

	for _, task := range localFailed {
		failFile(s.opts.Files.Tracker, task.trackerKey)
	}
	confirmed := taskRecords(localFailed)
	var remaining []*recordv1.Record
	if isRetryable(err) {
		remaining = pendingRecords
	} else {
		for _, task := range pending {
			failFile(s.opts.Files.Tracker, task.trackerKey)
		}
		confirmed = append(confirmed, pendingRecords...)
	}
	failures = append(failures, uploader.Failure{Message: err.Error(), Records: pendingRecords, Err: err})
	return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
}

// uploadLarge 走 prepare-multipart → 分片 PUT → complete-multipart 的大文件路径。
func (s *SaveSender) uploadLarge(ctx context.Context, tasks []*saveTask) uploader.Result {
	partSize := s.opts.Files.SavePart
	items := make([]api.PrepareMultipartSaveFile, len(tasks))
	for i, task := range tasks {
		items[i] = api.PrepareMultipartSaveFile{
			Path:     task.entry.path,
			Size:     task.entry.size,
			MD5:      task.entry.md5,
			MIMEType: task.entry.mimeType,
			Count:    int((task.entry.size + partSize - 1) / partSize),
		}
	}
	files, err := s.opts.API.PrepareMultipartSave(ctx, s.opts.ExperimentID, items)
	if err != nil {
		return batchFailure(taskRecords(tasks), err)
	}

	ready := make([]multipartReady, 0, len(tasks))
	storageFailed := make([]*saveTask, 0)
	localFailed := make([]*saveTask, 0)
	var failures []uploader.Failure
	for i, task := range tasks {
		trackFile(s.opts.Files.Tracker, task.trackerKey, task.entry.sourcePath, task.entry.size)
		parts, partErr := s.putParts(ctx, task, files[i].Parts, partSize)
		if partErr != nil {
			failures = append(failures, uploader.Failure{
				Message: saveMultipartFailedPrefix + task.entry.path,
				Records: []*recordv1.Record{task.record},
				Err:     partErr,
			})
			if putErrorRetryable(partErr) {
				storageFailed = append(storageFailed, task)
			} else {
				localFailed = append(localFailed, task)
			}
			continue
		}
		ready = append(ready, multipartReady{task: task, uploadID: files[i].UploadID, parts: parts})
	}

	if len(localFailed) > 0 {
		failedItems := make([]api.CompleteSaveFile, 0, len(localFailed))
		for _, task := range localFailed {
			failedItems = append(failedItems, api.CompleteSaveFile{Path: task.entry.path, State: saveStateFailed})
		}
		// 本地失败的上报失败不影响结果分类：记录已经终结。
		_ = s.opts.API.CompleteSaveFiles(ctx, s.opts.ExperimentID, failedItems)
		for _, task := range localFailed {
			failFile(s.opts.Files.Tracker, task.trackerKey)
		}
	}
	if len(ready) == 0 {
		return uploader.Result{
			Confirmed: taskRecords(localFailed),
			Remaining: taskRecords(storageFailed),
			Failures:  failures,
		}
	}

	completed := make([]api.CompletedMultipartSaveFile, len(ready))
	for i, item := range ready {
		completed[i] = api.CompletedMultipartSaveFile{
			Path:     item.task.entry.path,
			Parts:    item.parts,
			UploadID: item.uploadID,
		}
	}
	states, err := s.opts.API.CompleteMultipartSave(ctx, s.opts.ExperimentID, completed)
	if err != nil {
		return s.completeMultipartFailure(err, ready, storageFailed, localFailed, failures)
	}

	// 非空响应已由 API 层校验与请求一一对应且状态只含 UPLOADED/FAILED；
	// 空响应体沿用旧协议，视为全部成功。
	failedStates := make(map[string]struct{}, len(states))
	for _, state := range states {
		if state.State == api.MultipartFailed {
			failedStates[state.Path] = struct{}{}
		}
	}
	synced := 0
	confirmed := taskRecords(localFailed)
	for _, item := range ready {
		if _, failed := failedStates[item.task.entry.path]; failed {
			failures = append(failures, uploader.Failure{
				Message: saveMultipartStateFailedPrefix + item.task.entry.path,
				Records: []*recordv1.Record{item.task.record},
			})
			confirmed = append(confirmed, item.task.record)
			failFile(s.opts.Files.Tracker, item.task.trackerKey)
			continue
		}
		finishFile(s.opts.Files.Tracker, item.task.trackerKey)
		confirmed = append(confirmed, item.task.record)
		synced++
	}
	if synced > 0 {
		console.Infof("Synced %d save file(s)", synced)
	}
	return uploader.Result{Confirmed: confirmed, Remaining: taskRecords(storageFailed), Failures: failures}
}

// completeMultipartFailure 处理 complete-multipart 失败：可重试保留整批，否则终结。
func (s *SaveSender) completeMultipartFailure(
	err error,
	ready []multipartReady,
	storageFailed, localFailed []*saveTask,
	failures []uploader.Failure,
) uploader.Result {
	pending := make([]*saveTask, 0, len(ready)+len(storageFailed))
	for _, item := range ready {
		pending = append(pending, item.task)
	}
	pending = append(pending, storageFailed...)
	pendingRecords := taskRecords(pending)

	confirmed := taskRecords(localFailed)
	var remaining []*recordv1.Record
	if isRetryable(err) {
		remaining = pendingRecords
	} else {
		for _, task := range pending {
			failFile(s.opts.Files.Tracker, task.trackerKey)
		}
		confirmed = append(confirmed, pendingRecords...)
	}
	failures = append(failures, uploader.Failure{Message: err.Error(), Records: pendingRecords, Err: err})
	return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
}

// putParts 并发上传单个文件的分片，返回完成分片或首个错误。
func (s *SaveSender) putParts(ctx context.Context, task *saveTask, parts []api.MultipartPart, partSize int64) ([]api.CompletedPart, error) {
	file, err := os.Open(task.entry.sourcePath)
	if err != nil {
		return nil, err
	}
	defer func() { _ = file.Close() }()

	completed := make([]api.CompletedPart, len(parts))
	errs := runConcurrently(ctx, partSlots, len(parts), func(i int) error {
		part, err := s.putPart(ctx, file, task, parts[i], partSize)
		completed[i] = part
		return err
	})
	for _, err := range errs {
		if err != nil {
			return nil, err
		}
	}
	return completed, nil
}

// putPart 上传单个分片；失败时把该分片的进度回滚到 0。
func (s *SaveSender) putPart(
	ctx context.Context,
	file *os.File,
	task *saveTask,
	part api.MultipartPart,
	partSize int64,
) (api.CompletedPart, error) {
	offset := int64(part.PartNumber-1) * partSize
	length := min(partSize, max(task.entry.size-offset, 0))
	partKey := fmt.Sprintf("%s:part:%d", task.trackerKey, part.PartNumber)

	reader := io.Reader(io.NewSectionReader(file, offset, length))
	if tracker := s.opts.Files.Tracker; tracker != nil && length > 0 {
		var current int64
		reader = &progressReader{reader: reader, onRead: func(n int64) {
			current += n
			tracker.UpdateFileProgress(task.trackerKey, partKey, current)
		}}
	}
	etag, err := s.opts.Files.Putter.PutObject(ctx, part.URL, reader, length, savePartContentType)
	if err != nil {
		resetFileProgress(s.opts.Files.Tracker, task.trackerKey, partKey)
		return api.CompletedPart{}, err
	}
	return api.CompletedPart{PartNumber: part.PartNumber, ETag: etag}, nil
}

// putSave 上传单个文件，并在读取过程中上报字节进度；失败时回滚进度。
func (s *SaveSender) putSave(ctx context.Context, task *saveTask, url string) error {
	file, err := os.Open(task.entry.sourcePath)
	if err != nil {
		return err
	}
	defer func() { _ = file.Close() }()

	reader := io.Reader(file)
	if tracker := s.opts.Files.Tracker; tracker != nil && task.entry.size > 0 {
		var current int64
		reader = &progressReader{reader: reader, onRead: func(n int64) {
			current += n
			tracker.UpdateFileProgress(task.trackerKey, task.trackerKey, current)
		}}
	}
	if _, err := s.opts.Files.Putter.PutObject(ctx, url, reader, task.entry.size, task.entry.mimeType); err != nil {
		resetFileProgress(s.opts.Files.Tracker, task.trackerKey, task.trackerKey)
		return err
	}
	return nil
}

// resolveSource 解析 save 的本地文件：源路径优先，其次 run 目录下的 files 镜像。
func (s *SaveSender) resolveSource(save *savev1.SaveRecord) string {
	primary := save.GetSourcePath()
	if isRegularFile(primary) {
		return primary
	}
	if s.opts.Files.SkipStore {
		return ""
	}
	var fallback string
	if save.GetType() == savev1.SaveType_SAVE_TYPE_CUSTOM {
		parts := fs.WindowsParts(save.GetName())
		fallback = filepath.Join(append([]string{s.opts.Files.FilesDir}, parts...)...)
	} else {
		fallback = filepath.Join(s.opts.Files.FilesDir, fs.WindowsBase(save.GetSourcePath()))
	}
	if isRegularFile(fallback) {
		return fallback
	}
	return ""
}

// saveTrackerKey 构造 save 文件在 tracker 中的注册 key。
func saveTrackerKey(entry saveEntry) string {
	return fmt.Sprintf("%s:%d:%s", entry.path, entry.size, entry.md5)
}

// taskRecords 提取任务对应的记录。
func taskRecords(tasks []*saveTask) []*recordv1.Record {
	if len(tasks) == 0 {
		return nil
	}
	records := make([]*recordv1.Record, 0, len(tasks))
	for _, task := range tasks {
		records = append(records, task.record)
	}
	return records
}
