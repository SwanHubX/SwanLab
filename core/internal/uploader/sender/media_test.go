package sender

import (
	"context"
	"encoding/json"
	"io"
	"io/fs"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

const presignedPath = "/api/resources/presigned/put"

// objectPut 是一次假直传的记录。
type objectPut struct {
	url         string
	contentType string
	body        []byte
}

// fakePutter 在内存里接收 PUT，并按需返回失败。
type fakePutter struct {
	mu     sync.Mutex
	puts   []objectPut
	failOn func(url string) error
}

// PutObject 实现 ObjectPutter。
func (f *fakePutter) PutObject(_ context.Context, url string, body io.Reader, _ int64, contentType string) (string, error) {
	if f.failOn != nil {
		if err := f.failOn(url); err != nil {
			return "", err
		}
	}
	data, err := io.ReadAll(body)
	if err != nil {
		return "", err
	}
	f.mu.Lock()
	f.puts = append(f.puts, objectPut{url: url, contentType: contentType, body: data})
	f.mu.Unlock()
	return "etag", nil
}

// snapshot 返回已记录的 PUT。
func (f *fakePutter) snapshot() []objectPut {
	f.mu.Lock()
	defer f.mu.Unlock()
	return append([]objectPut(nil), f.puts...)
}

// withDefaults 填好 media sender 需要的目录与直传依赖。
func newMediaOptions(t *testing.T, cloud *clienttest.Cloud, files FileOptions) Options {
	t.Helper()
	opts := newSenderOptions(t, cloud)
	if files.MediaDir == "" {
		files.MediaDir = t.TempDir()
	}
	opts.Files = files
	return opts
}

// servePresigned 让假后端按请求里的 paths 返回等长 URL。
func servePresigned(cloud *clienttest.Cloud) {
	cloud.Handle(http.MethodPost, presignedPath, func(w http.ResponseWriter, _ *http.Request) {
		requests := cloud.Requests()
		var body struct {
			Paths []string `json:"paths"`
		}
		_ = json.Unmarshal(requests[len(requests)-1].Body, &body)
		urls := make([]string, len(body.Paths))
		for i, path := range body.Paths {
			urls[i] = "https://storage.test/" + path
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"urls": urls})
	})
}

// mediaRecord 构造一条媒体记录。
func mediaRecord(num int64, key string, step int64, typ columnv1.ColumnType, items ...*datav1.MediaItem) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Media{Media: &datav1.MediaRecord{
			Key:       key,
			Step:      step,
			Type:      typ,
			Timestamp: testTimestamp(),
			Value:     &datav1.MediaValue{Items: items},
		}},
	}
}

// writeMediaFile 在媒体目录下写一个本地文件。
func writeMediaFile(t *testing.T, dir, medium, name, content string) {
	t.Helper()
	target := filepath.Join(dir, medium, name)
	if err := os.MkdirAll(filepath.Dir(target), 0o755); err != nil {
		t.Fatalf("create media dir: %v", err)
	}
	if err := os.WriteFile(target, []byte(content), 0o644); err != nil {
		t.Fatalf("write media file: %v", err)
	}
}

// TestMediaUploadsPayloadAndFile 覆盖内存 payload 与本地文件两条收集路径，
// 以及指标字段、远端路径与内容类型。
func TestMediaUploadsPayloadAndFile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	servePresigned(cloud)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)

	mediaDir := t.TempDir()
	writeMediaFile(t, mediaDir, "image", "local.png", "local-bytes")
	putter := &fakePutter{}
	tracker := uploader.NewTracker()
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{
		MediaDir: mediaDir,
		Putter:   putter,
		Tracker:  tracker,
	}))

	record := mediaRecord(1, "train/sample", 7, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		&datav1.MediaItem{Filename: "payload.png", Payload: []byte("payload-bytes"), Caption: "first"},
		&datav1.MediaItem{Filename: "local.png"},
	)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}

	presignRequests := cloud.RequestsFor(http.MethodPost, presignedPath)
	if len(presignRequests) != 1 {
		t.Fatalf("presign requests = %d, want 1", len(presignRequests))
	}
	paths := decodeBody(t, presignRequests[0].Body)["paths"].([]any)
	wantPaths := []string{"media/image/payload.png", "media/image/local.png"}
	if len(paths) != len(wantPaths) {
		t.Fatalf("paths = %v, want %v", paths, wantPaths)
	}
	for i, want := range wantPaths {
		if paths[i] != want {
			t.Fatalf("paths[%d] = %v, want %q", i, paths[i], want)
		}
	}

	puts := putter.snapshot()
	if len(puts) != 2 {
		t.Fatalf("puts = %d, want 2", len(puts))
	}
	byURL := make(map[string]objectPut, len(puts))
	for _, put := range puts {
		byURL[put.url] = put
	}
	if got := byURL["https://storage.test/media/image/payload.png"]; string(got.body) != "payload-bytes" || got.contentType != "image/png" {
		t.Fatalf("payload put = %+v", got)
	}
	if got := byURL["https://storage.test/media/image/local.png"]; string(got.body) != "local-bytes" || got.contentType != "image/png" {
		t.Fatalf("file put = %+v", got)
	}

	metricRequests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(metricRequests) != 1 {
		t.Fatalf("metric requests = %d, want 1", len(metricRequests))
	}
	body := decodeBody(t, metricRequests[0].Body)
	if body["type"] != "media" {
		t.Fatalf("type = %v, want media", body["type"])
	}
	metrics := body["metrics"].([]any)
	if len(metrics) != 1 {
		t.Fatalf("metrics = %d, want 1", len(metrics))
	}
	metric := metrics[0].(map[string]any)
	if metric["key"] != "train/sample" || metric["index"] != float64(7) {
		t.Fatalf("metric = %+v", metric)
	}
	if metric["create_time"] != "2026-03-10T19:02:03.123400Z" {
		t.Fatalf("create_time = %v", metric["create_time"])
	}
	data := metric["data"].([]any)
	if len(data) != 2 || data[0] != wantPaths[0] || data[1] != wantPaths[1] {
		t.Fatalf("data = %v, want %v", data, wantPaths)
	}
	more := metric["more"].([]any)
	if len(more) != 1 || more[0].(map[string]any)["caption"] != "first" {
		t.Fatalf("more = %v, want one caption", more)
	}

	stats := tracker.Snapshot()
	if stats.TotalNumber != 2 || stats.UploadedNumber != 2 || len(stats.Files) != 0 {
		t.Fatalf("tracker = total %d uploaded %d files %d", stats.TotalNumber, stats.UploadedNumber, len(stats.Files))
	}
}

// TestMediaMissingLocalFile 覆盖本地缺失：确定性失败，且不申请预签名。
func TestMediaMissingLocalFile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	putter := &fakePutter{}
	tracker := uploader.NewTracker()
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: putter, Tracker: tracker}))

	record := mediaRecord(1, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		&datav1.MediaItem{Filename: "missing.png"},
	)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	if !strings.HasPrefix(result.Failures[0].Message, mediaFileNotFoundReason) {
		t.Fatalf("failure = %q, want %q prefix", result.Failures[0].Message, mediaFileNotFoundReason)
	}
	if got := cloud.CountFor(http.MethodPost, presignedPath); got != 0 {
		t.Fatalf("presign requests = %d, want 0", got)
	}
	if got := len(putter.snapshot()); got != 0 {
		t.Fatalf("puts = %d, want 0", got)
	}
	if stats := tracker.Snapshot(); stats.TotalNumber != 0 {
		t.Fatalf("tracker total = %d, want 0", stats.TotalNumber)
	}
}

// TestMediaSkipStoreSkipsLocalFile 覆盖 skip_store：即使本地有副本也不读取。
func TestMediaSkipStoreSkipsLocalFile(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	mediaDir := t.TempDir()
	writeMediaFile(t, mediaDir, "image", "present.png", "bytes")
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{
		MediaDir:  mediaDir,
		SkipStore: true,
		Putter:    &fakePutter{},
	}))

	record := mediaRecord(1, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		&datav1.MediaItem{Filename: "present.png"},
	)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}
	if !strings.Contains(result.Failures[0].Message, "skip_store") {
		t.Fatalf("failure = %q, want skip_store reason", result.Failures[0].Message)
	}
	if got := cloud.CountFor(http.MethodPost, presignedPath); got != 0 {
		t.Fatalf("presign requests = %d, want 0", got)
	}
}

// TestMediaPartialItemsUploadOnlyAvailable 覆盖一条记录里部分文件缺失：
// 只传可用文件，记录仍算成功。
func TestMediaPartialItemsUploadOnlyAvailable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	servePresigned(cloud)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	putter := &fakePutter{}
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: putter}))

	record := mediaRecord(1, "img", 2, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		&datav1.MediaItem{Filename: "ok.png", Payload: []byte("ok")},
		&datav1.MediaItem{Filename: "missing.png"},
	)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}

	paths := decodeBody(t, cloud.RequestsFor(http.MethodPost, presignedPath)[0].Body)["paths"].([]any)
	if len(paths) != 1 || paths[0] != "media/image/ok.png" {
		t.Fatalf("paths = %v, want only the available item", paths)
	}
	metrics := decodeBody(t, cloud.RequestsFor(http.MethodPost, metricsPath)[0].Body)["metrics"].([]any)
	data := metrics[0].(map[string]any)["data"].([]any)
	if len(data) != 1 || data[0] != "media/image/ok.png" {
		t.Fatalf("data = %v, want one path", data)
	}
	if got := len(putter.snapshot()); got != 1 {
		t.Fatalf("puts = %d, want 1", got)
	}
}

// TestMediaPutRetryableKeepsRecord 覆盖对象存储 403：过期链接保留重试，
// 下一轮会重新申请 URL。
func TestMediaPutRetryableKeepsRecord(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	servePresigned(cloud)
	tracker := uploader.NewTracker()
	putter := &fakePutter{failOn: func(string) error {
		return &client.APIError{Method: http.MethodPut, StatusCode: http.StatusForbidden}
	}}
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: putter, Tracker: tracker}))

	record := mediaRecord(1, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		&datav1.MediaItem{Filename: "a.png", Payload: []byte("x")},
	)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 0 || len(result.Remaining) != 1 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	if got := cloud.CountFor(http.MethodPost, metricsPath); got != 0 {
		t.Fatalf("metric requests = %d, want 0", got)
	}
	stats := tracker.Snapshot()
	if stats.TotalNumber != 1 || stats.UploadedNumber != 0 || len(stats.Files) != 1 {
		t.Fatalf("tracker = total %d uploaded %d files %d", stats.TotalNumber, stats.UploadedNumber, len(stats.Files))
	}
}

// TestMediaPutLocalErrorTerminates 覆盖本地读失败：重试不会恢复，终结该记录。
func TestMediaPutLocalErrorTerminates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	servePresigned(cloud)
	tracker := uploader.NewTracker()
	putter := &fakePutter{failOn: func(string) error {
		return &fs.PathError{Op: "open", Path: "a.png", Err: fs.ErrNotExist}
	}}
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: putter, Tracker: tracker}))

	record := mediaRecord(1, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		&datav1.MediaItem{Filename: "a.png", Payload: []byte("x")},
	)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	if got := cloud.CountFor(http.MethodPost, metricsPath); got != 0 {
		t.Fatalf("metric requests = %d, want 0", got)
	}
	stats := tracker.Snapshot()
	if stats.TotalNumber != 1 || stats.UploadedNumber != 0 || len(stats.Files) != 0 {
		t.Fatalf("tracker = total %d uploaded %d files %d", stats.TotalNumber, stats.UploadedNumber, len(stats.Files))
	}
}

// TestMediaMetricsFailureClassified 覆盖指标提交失败：5xx 保留重试，202 终结。
func TestMediaMetricsFailureClassified(t *testing.T) {
	cases := []struct {
		name      string
		status    int
		remaining int
		confirmed int
	}{
		{"retryable", http.StatusServiceUnavailable, 1, 0},
		{"house rejected", http.StatusAccepted, 0, 1},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			servePresigned(cloud)
			cloud.ServeJSON(http.MethodPost, metricsPath, tc.status, nil)
			tracker := uploader.NewTracker()
			sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: &fakePutter{}, Tracker: tracker}))

			record := mediaRecord(1, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
				&datav1.MediaItem{Filename: "a.png", Payload: []byte("x")},
			)
			result := sender.Upload(context.Background(), []*recordv1.Record{record})
			if len(result.Remaining) != tc.remaining || len(result.Confirmed) != tc.confirmed || len(result.Failures) != 1 {
				t.Fatalf("result = confirmed %d remaining %d failures %d",
					len(result.Confirmed), len(result.Remaining), len(result.Failures))
			}
			stats := tracker.Snapshot()
			if stats.TotalNumber != 1 || stats.UploadedNumber != 0 {
				t.Fatalf("tracker = total %d uploaded %d", stats.TotalNumber, stats.UploadedNumber)
			}
		})
	}
}

// TestMediaPresignFailureClassified 覆盖整批预签名失败的分类。
func TestMediaPresignFailureClassified(t *testing.T) {
	cases := []struct {
		name      string
		status    int
		remaining int
		confirmed int
	}{
		{"retryable", http.StatusServiceUnavailable, 1, 0},
		{"forbidden", http.StatusForbidden, 0, 1},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			cloud.ServeJSON(http.MethodPost, presignedPath, tc.status, map[string]string{"code": "nope"})
			tracker := uploader.NewTracker()
			sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: &fakePutter{}, Tracker: tracker}))

			record := mediaRecord(1, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
				&datav1.MediaItem{Filename: "a.png", Payload: []byte("x")},
			)
			result := sender.Upload(context.Background(), []*recordv1.Record{record})
			if len(result.Remaining) != tc.remaining || len(result.Confirmed) != tc.confirmed || len(result.Failures) != 1 {
				t.Fatalf("result = confirmed %d remaining %d failures %d",
					len(result.Confirmed), len(result.Remaining), len(result.Failures))
			}
			stats := tracker.Snapshot()
			wantFiles := 1
			if tc.remaining == 0 {
				wantFiles = 0
			}
			if len(stats.Files) != wantFiles {
				t.Fatalf("open files = %d, want %d", len(stats.Files), wantFiles)
			}
		})
	}
}

// TestMediaBlockedRecords 覆盖缺少 media 载荷与时间戳的确定性失败。
func TestMediaBlockedRecords(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: &fakePutter{}}))

	noPayload := &recordv1.Record{Num: 1}
	noTimestamp := mediaRecord(2, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		&datav1.MediaItem{Filename: "a.png", Payload: []byte("x")},
	)
	noTimestamp.GetMedia().Timestamp = nil

	result := sender.Upload(context.Background(), []*recordv1.Record{noPayload, noTimestamp})
	if len(result.Confirmed) != 2 || len(result.Remaining) != 0 || len(result.Failures) != 2 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	if result.Failures[0].Message != mediaMissingPayloadMessage || result.Failures[1].Message != mediaMissingTimestampMessage {
		t.Fatalf("failures = %+v", result.Failures)
	}
	if got := cloud.CountFor(http.MethodPost, presignedPath); got != 0 {
		t.Fatalf("presign requests = %d, want 0", got)
	}
}

// gatedPutter 把所有 PUT 挡在 release 上，用于观测并发峰值。
type gatedPutter struct {
	mu      sync.Mutex
	active  int
	peak    int
	entered chan struct{}
	release chan struct{}
}

// PutObject 实现 ObjectPutter：进入时记录峰值，随后阻塞到 release 关闭。
func (p *gatedPutter) PutObject(ctx context.Context, _ string, body io.Reader, _ int64, _ string) (string, error) {
	p.mu.Lock()
	p.active++
	p.peak = max(p.peak, p.active)
	p.mu.Unlock()

	select {
	case p.entered <- struct{}{}:
	default:
	}
	select {
	case <-p.release:
	case <-ctx.Done():
		return "", ctx.Err()
	}

	p.mu.Lock()
	p.active--
	p.mu.Unlock()
	_, err := io.Copy(io.Discard, body)
	return "", err
}

// maxActive 返回观测到的并发峰值。
func (p *gatedPutter) maxActive() int {
	p.mu.Lock()
	defer p.mu.Unlock()
	return p.peak
}

// TestMediaConcurrencyCap 覆盖进程级并发上限：20 个文件同时提交时峰值恰好为 10。
func TestMediaConcurrencyCap(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	servePresigned(cloud)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)

	putter := &gatedPutter{
		entered: make(chan struct{}, 20),
		release: make(chan struct{}),
	}
	sender := NewMedia(newMediaOptions(t, cloud, FileOptions{Putter: putter}))

	items := make([]*datav1.MediaItem, 0, 20)
	for i := 0; i < 20; i++ {
		items = append(items, &datav1.MediaItem{Filename: "f" + string(rune('a'+i)) + ".png", Payload: []byte("x")})
	}
	record := mediaRecord(1, "img", 1, columnv1.ColumnType_COLUMN_TYPE_IMAGE, items...)

	done := make(chan uploader.Result, 1)
	go func() { done <- sender.Upload(context.Background(), []*recordv1.Record{record}) }()

	for i := 0; i < mediaConcurrency; i++ {
		select {
		case <-putter.entered:
		case result := <-done:
			t.Fatalf("upload finished early: %+v", result)
		}
	}
	select {
	case <-putter.entered:
		t.Fatalf("more than %d PUTs entered concurrently", mediaConcurrency)
	case <-time.After(200 * time.Millisecond):
	}

	close(putter.release)
	result := <-done
	if len(result.Confirmed) != 1 || len(result.Failures) != 0 {
		t.Fatalf("result = %+v, want one confirmed record", result)
	}
	if got := putter.maxActive(); got != mediaConcurrency {
		t.Fatalf("peak concurrency = %d, want %d", got, mediaConcurrency)
	}
}
