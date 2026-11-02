package sender

import (
	"context"
	"crypto/md5"
	"encoding/hex"
	"errors"
	"io"
	"os"
	"sync"

	"golang.org/x/sync/semaphore"
	"google.golang.org/protobuf/encoding/protojson"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// md5Chunk 是计算 MD5 时的分块大小。
const md5Chunk = 8 << 20

// isRetryable 判断错误是否归 Remaining；父 ctx 取消时由 worker 停止调度。
func isRetryable(err error) bool {
	return errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) || client.IsRetryable(err)
}

// putErrorRetryable 判断对象存储 PUT 失败是否保留重试。
//
// 有 HTTP 响应（含 4xx 过期链接）或传输错误时，下一轮 prepare 会换新链接；
// 响应无法按协议解读或本地读写失败则终结；调用方取消保留给上层停止调度。
func putErrorRetryable(err error) bool {
	var apiErr *client.APIError
	if errors.As(err, &apiErr) {
		return !client.IsInvalidResponse(err)
	}
	return errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded)
}

// runConcurrently 在全局信号量下并发执行 fn，错误按索引写独立槽位。
func runConcurrently(ctx context.Context, slots *semaphore.Weighted, count int, fn func(int) error) []error {
	errs := make([]error, count)
	var group sync.WaitGroup
	for i := 0; i < count; i++ {
		if err := slots.Acquire(ctx, 1); err != nil {
			errs[i] = err
			continue
		}
		group.Add(1)
		go func() {
			defer group.Done()
			defer slots.Release(1)
			errs[i] = fn(i)
		}()
	}
	group.Wait()
	return errs
}

// progressReader 在每次读取后回调本次读取的字节数。
type progressReader struct {
	reader io.Reader
	onRead func(int64)
}

// Read 实现 io.Reader。
func (r *progressReader) Read(p []byte) (int, error) {
	n, err := r.reader.Read(p)
	if n > 0 && r.onRead != nil {
		r.onRead(int64(n))
	}
	return n, err
}

// withProgress 在需要时包装读取器，按累计字节上报文件或分片进度。
func withProgress(reader io.Reader, tracker *uploader.Tracker, fileKey, uploadKey string, size int64) io.Reader {
	if tracker == nil || size <= 0 {
		return reader
	}
	var current int64
	return &progressReader{reader: reader, onRead: func(n int64) {
		current += n
		tracker.UpdateFileProgress(fileKey, uploadKey, current)
	}}
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

// computeMD5s 并发计算一批文件的 MD5，结果与错误按下标对齐。
func computeMD5s(ctx context.Context, paths []string) ([]string, []error) {
	digests := make([]string, len(paths))
	errs := runConcurrently(ctx, md5Slots, len(paths), func(i int) error {
		digest, err := fileMD5(paths[i])
		digests[i] = digest
		return err
	})
	return digests, errs
}

// isRegularFile 判断路径是普通文件；不存在或不可访问时返回 false。
func isRegularFile(path string) bool {
	if path == "" {
		return false
	}
	info, err := os.Stat(path)
	return err == nil && !info.IsDir()
}

// encodeTimestamp 复用 protojson 的 0/3/6/9 位小数排版，去掉外层引号。
func encodeTimestamp(ts *timestamppb.Timestamp) (string, bool) {
	raw, err := protojson.Marshal(ts)
	if err != nil || len(raw) < 2 || raw[0] != '"' {
		return "", false
	}
	return string(raw[1 : len(raw)-1]), true
}

// trackFile 注册文件进度行；nil tracker 与 size<=0 都忽略。
func trackFile(tracker *uploader.Tracker, key, path string, size int64) {
	if tracker == nil || size <= 0 {
		return
	}
	tracker.AddFile(key, path, size)
}

// finishFile 标记文件上传完成。
func finishFile(tracker *uploader.Tracker, key string) {
	if tracker != nil {
		tracker.FinishFile(key)
	}
}

// failFile 关闭未成功文件的进度行。
func failFile(tracker *uploader.Tracker, key string) {
	if tracker != nil {
		tracker.FailFile(key)
	}
}

// resetFileProgress 把文件或分片已上报的进度回滚到 0，失败重试时避免重复计字节。
func resetFileProgress(tracker *uploader.Tracker, fileKey, uploadKey string) {
	if tracker != nil {
		tracker.UpdateFileProgress(fileKey, uploadKey, 0)
	}
}

// blockedFailure 包装确定性失败的记录。
func blockedFailure(message string, records ...*recordv1.Record) *uploader.Failure {
	return &uploader.Failure{Message: message, Records: records}
}

// blockedFailures 包装一批确定性失败的记录；空集合返回 nil。
func blockedFailures(message string, records []*recordv1.Record) []uploader.Failure {
	if len(records) == 0 {
		return nil
	}
	return []uploader.Failure{*blockedFailure(message, records...)}
}

// confirmedFailure 构造「记录已终结但带失败诊断」的结果。
func confirmedFailure(record *recordv1.Record, message string) uploader.Result {
	records := []*recordv1.Record{record}
	return uploader.Result{Confirmed: records, Failures: []uploader.Failure{*blockedFailure(message, record)}}
}

// batchFailure 把整批请求失败按可重试语义分流。
func batchFailure(records []*recordv1.Record, err error) uploader.Result {
	failure := uploader.Failure{Message: err.Error(), Records: records, Err: err}
	if isRetryable(err) {
		return uploader.Result{Remaining: records, Failures: []uploader.Failure{failure}}
	}
	return uploader.Result{Confirmed: records, Failures: []uploader.Failure{failure}}
}

// mergeResult 把一次子上传的结果并入累计结果。
func mergeResult(base, next uploader.Result) uploader.Result {
	base.Confirmed = append(base.Confirmed, next.Confirmed...)
	base.Remaining = append(base.Remaining, next.Remaining...)
	base.Failures = append(base.Failures, next.Failures...)
	return base
}
