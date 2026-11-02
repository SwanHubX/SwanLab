package sender

import (
	"context"
	"crypto/md5"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"unicode/utf8"

	"golang.org/x/sync/errgroup"
	"gopkg.in/yaml.v3"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
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
	// saveInvalidBatchPrefix 是非法分批配置的终结原因前缀。
	saveInvalidBatchPrefix = "invalid save batch size "
	// saveInternalUnavailablePrefix 是内部 save 内容不可用的终结原因前缀。
	saveInternalUnavailablePrefix = "internal save content unavailable: "
	// saveInternalUnparseablePrefix 是内部 save 解析失败的终结原因前缀。
	saveInternalUnparseablePrefix = "internal save payload unparseable: "

	// saveStateUploaded 是 complete 阶段上报的成功状态。
	saveStateUploaded = "UPLOADED"
	// saveStateFailed 是 complete 阶段上报的本地失败状态。
	saveStateFailed = "FAILED"
	// md5Chunk 是计算 MD5 时的分块大小。
	md5Chunk = 8 << 20
)

// SaveSender 分发内部元文件与 CUSTOM 文件：前者进 Profile，后者直传对象存储。
type SaveSender struct{ opts Options }

// NewSave 构造 save sender。
func NewSave(opts Options) *SaveSender { return &SaveSender{opts: opts} }

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
		var content map[string]any
		if err := json.Unmarshal(data, &content); err != nil || content == nil {
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
	if batch <= 0 {
		batchRecords := pendingRecords(pending)
		failures = append(failures, uploader.Failure{
			Message: fmt.Sprintf("%s%d", saveInvalidBatchPrefix, batch),
			Records: batchRecords,
		})
		return uploader.Result{Confirmed: append(confirmed, batchRecords...), Failures: failures}
	}

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
	digests, errs := computeMD5s(ctx, pending)
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
	result := s.uploadSmall(ctx, tasks)
	result.Confirmed = append(confirmed, result.Confirmed...)
	result.Remaining = append(remaining, result.Remaining...)
	result.Failures = append(failures, result.Failures...)
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
		s.trackFile(task)
	}

	putErrors := s.putSaves(ctx, tasks, urls)
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
		s.finishFile(task)
	}
	for _, task := range localFailed {
		s.failFile(task)
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
		s.failFile(task)
	}
	confirmed := taskRecords(localFailed)
	var remaining []*recordv1.Record
	if isRetryable(err) {
		remaining = pendingRecords
	} else {
		for _, task := range pending {
			s.failFile(task)
		}
		confirmed = append(confirmed, pendingRecords...)
	}
	failures = append(failures, uploader.Failure{Message: err.Error(), Records: pendingRecords, Err: err})
	return uploader.Result{Confirmed: confirmed, Remaining: remaining, Failures: failures}
}

// putSaves 并发直传小文件，错误按索引写独立槽位。
func (s *SaveSender) putSaves(ctx context.Context, tasks []*saveTask, urls []string) []error {
	errs := make([]error, len(tasks))
	var group errgroup.Group
	for i := range tasks {
		if err := saveSlots.Acquire(ctx, 1); err != nil {
			errs[i] = err
			continue
		}
		group.Go(func() error {
			defer saveSlots.Release(1)
			errs[i] = s.putSave(ctx, tasks[i], urls[i])
			return nil
		})
	}
	_ = group.Wait()
	return errs
}

// putSave 上传单个文件，并在读取过程中上报字节进度。
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
	_, err = s.opts.Files.Putter.PutObject(ctx, url, reader, task.entry.size, task.entry.mimeType)
	return err
}

// computeMD5s 并发计算 MD5，结果与错误按下标对齐。
func computeMD5s(ctx context.Context, pending []pendingSave) ([]string, []error) {
	digests := make([]string, len(pending))
	errs := make([]error, len(pending))
	var group errgroup.Group
	for i := range pending {
		if err := md5Slots.Acquire(ctx, 1); err != nil {
			errs[i] = err
			continue
		}
		group.Go(func() error {
			defer md5Slots.Release(1)
			digests[i], errs[i] = fileMD5(pending[i].sourcePath)
			return nil
		})
	}
	_ = group.Wait()
	return digests, errs
}

// fileMD5 分块计算文件 MD5，避免大文件整块读入内存。
func fileMD5(path string) (string, error) {
	file, err := os.Open(path)
	if err != nil {
		return "", err
	}
	defer func() { _ = file.Close() }()

	hash := md5.New()
	if _, err := io.CopyBuffer(hash, file, make([]byte, md5Chunk)); err != nil {
		return "", err
	}
	return hex.EncodeToString(hash.Sum(nil)), nil
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
		parts := windowsParts(save.GetName())
		fallback = filepath.Join(append([]string{s.opts.Files.FilesDir}, parts...)...)
	} else {
		fallback = filepath.Join(s.opts.Files.FilesDir, windowsBase(save.GetSourcePath()))
	}
	if isRegularFile(fallback) {
		return fallback
	}
	return ""
}

// isRegularFile 判断路径是普通文件；不存在或不可访问时返回 false。
func isRegularFile(path string) bool {
	if path == "" {
		return false
	}
	info, err := os.Stat(path)
	return err == nil && !info.IsDir()
}

// windowsParts 按 Windows 语义拆分相对路径："/" 与 "\" 都是分隔符。
func windowsParts(name string) []string {
	parts := strings.FieldsFunc(name, func(r rune) bool { return r == '/' || r == '\\' })
	kept := parts[:0]
	for _, part := range parts {
		if part == "." {
			continue
		}
		kept = append(kept, part)
	}
	return kept
}

// windowsBase 返回 Windows 语义下的最后一段文件名。
func windowsBase(path string) string {
	if index := strings.LastIndexAny(path, `/\`); index >= 0 {
		return path[index+1:]
	}
	return path
}

// saveTrackerKey 构造 save 文件在 tracker 中的注册 key。
func saveTrackerKey(entry saveEntry) string {
	return fmt.Sprintf("%s:%d:%s", entry.path, entry.size, entry.md5)
}

// trackFile 在上传前注册进度行；size<=0 由 tracker 自行忽略。
func (s *SaveSender) trackFile(task *saveTask) {
	if s.opts.Files.Tracker == nil {
		return
	}
	s.opts.Files.Tracker.AddFile(task.trackerKey, task.entry.sourcePath, task.entry.size)
}

// finishFile 标记文件上传完成。
func (s *SaveSender) finishFile(task *saveTask) {
	if s.opts.Files.Tracker != nil {
		s.opts.Files.Tracker.FinishFile(task.trackerKey)
	}
}

// failFile 关闭未成功文件的进度行。
func (s *SaveSender) failFile(task *saveTask) {
	if s.opts.Files.Tracker != nil {
		s.opts.Files.Tracker.FailFile(task.trackerKey)
	}
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

// pendingRecords 提取待上传记录。
func pendingRecords(pending []pendingSave) []*recordv1.Record {
	records := make([]*recordv1.Record, 0, len(pending))
	for _, item := range pending {
		records = append(records, item.record)
	}
	return records
}

// batchFailure 把整批请求失败按可重试语义分流。
func batchFailure(records []*recordv1.Record, err error) uploader.Result {
	failure := uploader.Failure{Message: err.Error(), Records: records, Err: err}
	if isRetryable(err) {
		return uploader.Result{Remaining: records, Failures: []uploader.Failure{failure}}
	}
	return uploader.Result{Confirmed: records, Failures: []uploader.Failure{failure}}
}

// confirmedFailure 构造「记录已终结但带失败诊断」的结果。
func confirmedFailure(record *recordv1.Record, message string) uploader.Result {
	records := []*recordv1.Record{record}
	return uploader.Result{Confirmed: records, Failures: []uploader.Failure{{Message: message, Records: records}}}
}
