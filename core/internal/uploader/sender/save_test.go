package sender

import (
	"context"
	"encoding/json"
	"fmt"
	"io"
	"io/fs"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

const (
	profilePath  = "/api/project/alice/demo/runs/exp-1/profile"
	preparePath  = "/api/experiment/exp-1/files/prepare"
	completePath = "/api/experiment/exp-1/files/complete"
)

// defaultSaveFiles 填好 CUSTOM 上传需要的目录、上限与分批。
func defaultSaveFiles(t *testing.T, putter ObjectPutter) FileOptions {
	t.Helper()
	return FileOptions{
		FilesDir:  t.TempDir(),
		SaveSize:  50 << 30,
		SaveSplit: 100 << 20,
		SavePart:  32 << 20,
		SaveBatch: 100,
		Putter:    putter,
	}
}

// newSaveSender 构造 save sender。
func newSaveSender(t *testing.T, cloud *clienttest.Cloud, files FileOptions) *SaveSender {
	t.Helper()
	opts := newSenderOptions(t, cloud)
	opts.Files = files
	return NewSave(opts)
}

// saveRecord 构造一条 save 记录；payload 为 nil 表示字段缺失。
func saveRecord(num int64, typ savev1.SaveType, name, sourcePath string, payload []byte) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Save{Save: &savev1.SaveRecord{
			Name:       name,
			SourcePath: sourcePath,
			Type:       typ,
			Payload:    payload,
		}},
	}
}

// writeSaveFile 写一个本地文件并返回路径。
func writeSaveFile(t *testing.T, dir, name string, content []byte) string {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatalf("create dir: %v", err)
	}
	if err := os.WriteFile(path, content, 0o644); err != nil {
		t.Fatalf("write file: %v", err)
	}
	return path
}

// serveSavePrepare 让假后端按请求里的 files 返回等长 URL。
func serveSavePrepare(cloud *clienttest.Cloud) {
	cloud.Handle(http.MethodPost, preparePath, func(w http.ResponseWriter, _ *http.Request) {
		requests := cloud.Requests()
		var body struct {
			Files []map[string]any `json:"files"`
		}
		_ = json.Unmarshal(requests[len(requests)-1].Body, &body)
		urls := make([]string, len(body.Files))
		for i := range body.Files {
			urls[i] = fmt.Sprintf("https://storage.test/save/%d", i)
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"urls": urls})
	})
}

// serveSaveComplete 让 complete 返回成功。
func serveSaveComplete(cloud *clienttest.Cloud) {
	cloud.ServeJSON(http.MethodPost, completePath, http.StatusCreated, map[string]any{"ok": true})
}

// prepareBatches 返回每次 prepare 请求里的 files 列表。
func prepareBatches(t *testing.T, cloud *clienttest.Cloud) [][]map[string]any {
	t.Helper()
	requests := cloud.RequestsFor(http.MethodPost, preparePath)
	batches := make([][]map[string]any, 0, len(requests))
	for _, request := range requests {
		files, _ := decodeBody(t, request.Body)["files"].([]any)
		batch := make([]map[string]any, 0, len(files))
		for _, item := range files {
			batch = append(batch, item.(map[string]any))
		}
		batches = append(batches, batch)
	}
	return batches
}

// completeBatches 返回每次 complete 请求里的 files 列表。
func completeBatches(t *testing.T, cloud *clienttest.Cloud) [][]map[string]any {
	t.Helper()
	requests := cloud.RequestsFor(http.MethodPost, completePath)
	batches := make([][]map[string]any, 0, len(requests))
	for _, request := range requests {
		files, _ := decodeBody(t, request.Body)["files"].([]any)
		batch := make([]map[string]any, 0, len(files))
		for _, item := range files {
			batch = append(batch, item.(map[string]any))
		}
		batches = append(batches, batch)
	}
	return batches
}

// TestSaveInternalUploadsFourKinds 覆盖四类内部元文件解码与 Profile 上报。
func TestSaveInternalUploadsFourKinds(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPut, profilePath, http.StatusOK, nil)
	sender := newSaveSender(t, cloud, defaultSaveFiles(t, &fakePutter{}))

	records := []*recordv1.Record{
		saveRecord(1, savev1.SaveType_SAVE_TYPE_METADATA, "swanlab-metadata.json", "/missing/meta.json", []byte(`{"loss": 1}`)),
		saveRecord(2, savev1.SaveType_SAVE_TYPE_CONFIG, "config.yaml", "/missing/config.yaml", []byte("epochs: 3\n")),
		saveRecord(3, savev1.SaveType_SAVE_TYPE_REQUIREMENTS, "requirements.txt", "/missing/requirements.txt", []byte("numpy==1.0\n")),
		saveRecord(4, savev1.SaveType_SAVE_TYPE_CONDA, "conda.yaml", "/missing/conda.yaml", []byte("name: env\n")),
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 4 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}

	requests := cloud.RequestsFor(http.MethodPut, profilePath)
	if len(requests) != 4 {
		t.Fatalf("profile requests = %d, want 4", len(requests))
	}
	metadata, _ := decodeBody(t, requests[0].Body)["metadata"].(map[string]any)
	if metadata["loss"] != float64(1) {
		t.Fatalf("metadata = %+v", metadata)
	}
	config, _ := decodeBody(t, requests[1].Body)["config"].(map[string]any)
	if config["epochs"] != float64(3) {
		t.Fatalf("config = %+v", config)
	}
	if got := decodeBody(t, requests[2].Body)["requirements"]; got != "numpy==1.0\n" {
		t.Fatalf("requirements = %q", got)
	}
	if got := decodeBody(t, requests[3].Body)["conda"]; got != "name: env\n" {
		t.Fatalf("conda = %q", got)
	}
	if got := cloud.CountFor(http.MethodPost, preparePath); got != 0 {
		t.Fatalf("prepare requests = %d, want 0", got)
	}
}

// TestSaveInternalReadsFallbackSources 覆盖 payload 缺失时的源路径与 files 目录回退。
func TestSaveInternalReadsFallbackSources(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPut, profilePath, http.StatusOK, nil)
	files := defaultSaveFiles(t, &fakePutter{})
	sender := newSaveSender(t, cloud, files)

	primary := writeSaveFile(t, t.TempDir(), "meta.json", []byte(`{"source": "primary"}`))
	writeSaveFile(t, files.FilesDir, "config.yaml", []byte("seed: 7\n"))
	records := []*recordv1.Record{
		saveRecord(1, savev1.SaveType_SAVE_TYPE_METADATA, "swanlab-metadata.json", primary, nil),
		saveRecord(2, savev1.SaveType_SAVE_TYPE_CONFIG, "config.yaml", "/missing/config.yaml", nil),
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 2 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}

	requests := cloud.RequestsFor(http.MethodPut, profilePath)
	if len(requests) != 2 {
		t.Fatalf("profile requests = %d, want 2", len(requests))
	}
	metadata, _ := decodeBody(t, requests[0].Body)["metadata"].(map[string]any)
	if metadata["source"] != "primary" {
		t.Fatalf("metadata = %+v, want source from source_path", metadata)
	}
	config, _ := decodeBody(t, requests[1].Body)["config"].(map[string]any)
	if config["seed"] != float64(7) {
		t.Fatalf("config = %+v, want seed from files_dir", config)
	}
}

// TestSaveInternalSkipStoreDoesNotFallback 覆盖 skip_store：禁止回退到旧 run 的同名文件。
func TestSaveInternalSkipStoreDoesNotFallback(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	files := defaultSaveFiles(t, &fakePutter{})
	files.SkipStore = true
	writeSaveFile(t, files.FilesDir, "config.yaml", []byte("seed: 7\n"))
	sender := newSaveSender(t, cloud, files)

	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CONFIG, "config.yaml", "/missing/config.yaml", nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	if !strings.HasPrefix(result.Failures[0].Message, saveInternalUnavailablePrefix) {
		t.Fatalf("failure = %q, want %q prefix", result.Failures[0].Message, saveInternalUnavailablePrefix)
	}
	if got := cloud.CountFor(http.MethodPut, profilePath); got != 0 {
		t.Fatalf("profile requests = %d, want 0", got)
	}
}

// TestSaveInternalUnparseableTerminates 覆盖内部 save 解析失败：确定性终结且不发请求。
func TestSaveInternalUnparseableTerminates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	sender := newSaveSender(t, cloud, defaultSaveFiles(t, &fakePutter{}))

	records := []*recordv1.Record{
		saveRecord(1, savev1.SaveType_SAVE_TYPE_METADATA, "meta.json", "/missing/meta.json", []byte("not-json")),
		saveRecord(2, savev1.SaveType_SAVE_TYPE_CONFIG, "config.yaml", "/missing/config.yaml", []byte("- 1\n- 2\n")),
		saveRecord(3, savev1.SaveType_SAVE_TYPE_REQUIREMENTS, "requirements.txt", "/missing/requirements.txt", []byte{}),
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 3 || len(result.Remaining) != 0 || len(result.Failures) != 3 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	for _, failure := range result.Failures {
		if !strings.HasPrefix(failure.Message, saveInternalUnparseablePrefix) {
			t.Fatalf("failure = %q, want %q prefix", failure.Message, saveInternalUnparseablePrefix)
		}
	}
	if got := cloud.CountFor(http.MethodPut, profilePath); got != 0 {
		t.Fatalf("profile requests = %d, want 0", got)
	}
}

// TestSaveDispatchErrors 覆盖缺 save 载荷与 CUSTOM 携带 payload 的协议错误。
func TestSaveDispatchErrors(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	files := defaultSaveFiles(t, &fakePutter{})
	sender := newSaveSender(t, cloud, files)

	noPayload := &recordv1.Record{Num: 1}
	customWithPayload := saveRecord(2, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", "/missing/a.txt", []byte("x"))
	missing := saveRecord(3, savev1.SaveType_SAVE_TYPE_CUSTOM, "b.txt", "/missing/b.txt", nil)

	result := sender.Upload(context.Background(), []*recordv1.Record{noPayload, customWithPayload, missing})
	if len(result.Confirmed) != 3 || len(result.Remaining) != 0 || len(result.Failures) != 3 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	messages := []string{result.Failures[0].Message, result.Failures[1].Message, result.Failures[2].Message}
	want := []string{saveNoPayloadMessage, saveCustomPayloadMessage, saveFileNotFoundPrefix + "/missing/b.txt"}
	for i, prefix := range want {
		if messages[i] != prefix {
			t.Fatalf("failure[%d] = %q, want %q", i, messages[i], prefix)
		}
	}
	if got := cloud.CountFor(http.MethodPost, preparePath); got != 0 {
		t.Fatalf("prepare requests = %d, want 0", got)
	}
}

// TestSaveInternalRetryable 覆盖 Profile 上报失败分类。
func TestSaveInternalRetryable(t *testing.T) {
	cases := []struct {
		name      string
		status    int
		remaining int
		confirmed int
	}{
		{"retryable", http.StatusInternalServerError, 1, 0},
		{"bad request", http.StatusBadRequest, 0, 1},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			cloud.ServeJSON(http.MethodPut, profilePath, tc.status, map[string]string{"code": "nope"})
			sender := newSaveSender(t, cloud, defaultSaveFiles(t, &fakePutter{}))

			record := saveRecord(1, savev1.SaveType_SAVE_TYPE_METADATA, "meta.json", "/missing/meta.json", []byte(`{"a": 1}`))
			result := sender.Upload(context.Background(), []*recordv1.Record{record})
			if len(result.Remaining) != tc.remaining || len(result.Confirmed) != tc.confirmed || len(result.Failures) != 1 {
				t.Fatalf("result = confirmed %d remaining %d failures %d",
					len(result.Confirmed), len(result.Remaining), len(result.Failures))
			}
		})
	}
}

// TestSaveSmallUpload 覆盖 prepare → PUT → complete 的完整字段与 tracker。
func TestSaveSmallUpload(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveSavePrepare(cloud)
	serveSaveComplete(cloud)
	putter := &fakePutter{}
	tracker := uploader.NewTracker()
	files := defaultSaveFiles(t, putter)
	files.Tracker = tracker
	sender := newSaveSender(t, cloud, files)

	source := writeSaveFile(t, t.TempDir(), "artifact.py", []byte("hello world"))
	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "checkpoints/artifact.py", source, nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}

	batches := prepareBatches(t, cloud)
	if len(batches) != 1 || len(batches[0]) != 1 {
		t.Fatalf("prepare batches = %+v, want one file", batches)
	}
	item := batches[0][0]
	if item["path"] != "checkpoints/artifact.py" || item["size"] != float64(11) {
		t.Fatalf("prepare item = %+v", item)
	}
	if item["md5"] != "5eb63bbbe01eeed093cb22bb8f5acdc3" {
		t.Fatalf("md5 = %v, want md5 of hello world", item["md5"])
	}
	if item["mimeType"] != "text/x-python" {
		t.Fatalf("mimeType = %v, want text/x-python", item["mimeType"])
	}

	puts := putter.snapshot()
	if len(puts) != 1 {
		t.Fatalf("puts = %d, want 1", len(puts))
	}
	if string(puts[0].body) != "hello world" || puts[0].contentType != "text/x-python" {
		t.Fatalf("put = %+v", puts[0])
	}

	completed := completeBatches(t, cloud)
	if len(completed) != 1 || len(completed[0]) != 1 {
		t.Fatalf("complete batches = %+v", completed)
	}
	if completed[0][0]["path"] != "checkpoints/artifact.py" || completed[0][0]["state"] != saveStateUploaded {
		t.Fatalf("complete item = %+v", completed[0][0])
	}

	stats := tracker.Snapshot()
	if stats.TotalNumber != 1 || stats.UploadedNumber != 1 {
		t.Fatalf("tracker = total %d uploaded %d", stats.TotalNumber, stats.UploadedNumber)
	}
}

// TestSaveSmallBatches 覆盖按 SaveBatch 分批：105 个文件拆成 100 + 5。
func TestSaveSmallBatches(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveSavePrepare(cloud)
	serveSaveComplete(cloud)
	putter := &fakePutter{}
	sender := newSaveSender(t, cloud, defaultSaveFiles(t, putter))

	dir := t.TempDir()
	records := make([]*recordv1.Record, 0, 105)
	for i := 0; i < 105; i++ {
		name := fmt.Sprintf("f%03d.txt", i)
		source := writeSaveFile(t, dir, name, []byte("x"))
		records = append(records, saveRecord(int64(i+1), savev1.SaveType_SAVE_TYPE_CUSTOM, name, source, nil))
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 105 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}

	batches := prepareBatches(t, cloud)
	if len(batches) != 2 || len(batches[0]) != 100 || len(batches[1]) != 5 {
		t.Fatalf("prepare batch sizes = %v, want 100 + 5", batchSizes(batches))
	}
	if batches[0][0]["path"] != "f000.txt" || batches[0][99]["path"] != "f099.txt" {
		t.Fatalf("first batch order = %v .. %v", batches[0][0]["path"], batches[0][99]["path"])
	}
	if batches[1][0]["path"] != "f100.txt" || batches[1][4]["path"] != "f104.txt" {
		t.Fatalf("second batch order = %v .. %v", batches[1][0]["path"], batches[1][4]["path"])
	}
	if got := len(putter.snapshot()); got != 105 {
		t.Fatalf("puts = %d, want 105", got)
	}
	if got := len(completeBatches(t, cloud)); got != 2 {
		t.Fatalf("complete batches = %d, want 2", got)
	}
}

// batchSizes 返回每批的文件数。
func batchSizes(batches [][]map[string]any) []int {
	sizes := make([]int, 0, len(batches))
	for _, batch := range batches {
		sizes = append(sizes, len(batch))
	}
	return sizes
}

// TestFileMD5 覆盖已知向量与空文件。
func TestFileMD5(t *testing.T) {
	cases := []struct {
		name    string
		content []byte
		want    string
	}{
		{"abc", []byte("abc"), "900150983cd24fb0d6963f7d28e17f72"},
		{"empty", nil, "d41d8cd98f00b204e9800998ecf8427e"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			path := writeSaveFile(t, t.TempDir(), "file.bin", tc.content)
			got, err := fileMD5(path)
			if err != nil {
				t.Fatalf("fileMD5: %v", err)
			}
			if got != tc.want {
				t.Fatalf("md5 = %q, want %q", got, tc.want)
			}
		})
	}
	if _, err := fileMD5(filepath.Join(t.TempDir(), "missing")); err == nil {
		t.Fatal("fileMD5 on a missing file should fail")
	}
}

// TestSavePrepareFailure 覆盖 prepare 失败：4xx 终结、5xx 保留重试。
func TestSavePrepareFailure(t *testing.T) {
	cases := []struct {
		name      string
		status    int
		remaining int
		confirmed int
	}{
		{"forbidden", http.StatusForbidden, 0, 1},
		{"unavailable", http.StatusServiceUnavailable, 1, 0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			cloud.ServeJSON(http.MethodPost, preparePath, tc.status, map[string]string{"code": "nope"})
			putter := &fakePutter{}
			tracker := uploader.NewTracker()
			files := defaultSaveFiles(t, putter)
			files.Tracker = tracker
			sender := newSaveSender(t, cloud, files)

			source := writeSaveFile(t, t.TempDir(), "a.txt", []byte("x"))
			record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", source, nil)
			result := sender.Upload(context.Background(), []*recordv1.Record{record})
			if len(result.Confirmed) != tc.confirmed || len(result.Remaining) != tc.remaining || len(result.Failures) != 1 {
				t.Fatalf("result = confirmed %d remaining %d failures %d",
					len(result.Confirmed), len(result.Remaining), len(result.Failures))
			}
			if got := len(putter.snapshot()); got != 0 {
				t.Fatalf("puts = %d, want 0", got)
			}
			if stats := tracker.Snapshot(); stats.TotalNumber != 0 {
				t.Fatalf("tracker total = %d, want 0 (tracking starts after prepare)", stats.TotalNumber)
			}
		})
	}
}

// TestSavePutFailureClassification 覆盖单文件直传失败：对象存储错误保留重试，
// 本地错误终结，两者都上报 complete。
func TestSavePutFailureClassification(t *testing.T) {
	t.Run("object storage", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		serveSavePrepare(cloud)
		serveSaveComplete(cloud)
		putter := &fakePutter{failOn: func(url string) error {
			if strings.HasSuffix(url, "/save/0") {
				return &client.APIError{Method: http.MethodPut, StatusCode: http.StatusForbidden}
			}
			return nil
		}}
		sender := newSaveSender(t, cloud, defaultSaveFiles(t, putter))

		dir := t.TempDir()
		records := []*recordv1.Record{
			saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", writeSaveFile(t, dir, "a.txt", []byte("a")), nil),
			saveRecord(2, savev1.SaveType_SAVE_TYPE_CUSTOM, "b.txt", writeSaveFile(t, dir, "b.txt", []byte("b")), nil),
		}
		result := sender.Upload(context.Background(), records)
		if len(result.Confirmed) != 1 || len(result.Remaining) != 1 || len(result.Failures) != 1 {
			t.Fatalf("result = confirmed %d remaining %d failures %d",
				len(result.Confirmed), len(result.Remaining), len(result.Failures))
		}
		completed := completeBatches(t, cloud)
		if len(completed) != 1 || len(completed[0]) != 1 || completed[0][0]["path"] != "b.txt" {
			t.Fatalf("complete = %+v, want only the uploaded file", completed)
		}
	})

	t.Run("local read error", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		serveSavePrepare(cloud)
		serveSaveComplete(cloud)
		putter := &fakePutter{failOn: func(url string) error {
			if strings.HasSuffix(url, "/save/0") {
				return &fs.PathError{Op: "open", Path: "a.txt", Err: fs.ErrNotExist}
			}
			return nil
		}}
		tracker := uploader.NewTracker()
		files := defaultSaveFiles(t, putter)
		files.Tracker = tracker
		sender := newSaveSender(t, cloud, files)

		dir := t.TempDir()
		records := []*recordv1.Record{
			saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", writeSaveFile(t, dir, "a.txt", []byte("a")), nil),
			saveRecord(2, savev1.SaveType_SAVE_TYPE_CUSTOM, "b.txt", writeSaveFile(t, dir, "b.txt", []byte("b")), nil),
		}
		result := sender.Upload(context.Background(), records)
		if len(result.Confirmed) != 2 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
			t.Fatalf("result = confirmed %d remaining %d failures %d",
				len(result.Confirmed), len(result.Remaining), len(result.Failures))
		}
		completed := completeBatches(t, cloud)
		if len(completed) != 1 || len(completed[0]) != 2 {
			t.Fatalf("complete = %+v, want uploaded + failed", completed)
		}
		states := map[string]any{}
		for _, item := range completed[0] {
			states[item["path"].(string)] = item["state"]
		}
		if states["a.txt"] != saveStateFailed || states["b.txt"] != saveStateUploaded {
			t.Fatalf("complete states = %+v", states)
		}
		if stats := tracker.Snapshot(); stats.TotalNumber != 2 || stats.UploadedNumber != 1 || len(stats.Files) != 0 {
			t.Fatalf("tracker = total %d uploaded %d files %d", stats.TotalNumber, stats.UploadedNumber, len(stats.Files))
		}
	})
}

// TestSaveCompleteFailureKeepsStorageFailures 覆盖 complete 可重试时，
// 已上传与对象存储失败的记录都要保留。
func TestSaveCompleteFailureKeepsStorageFailures(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveSavePrepare(cloud)
	cloud.ServeJSON(http.MethodPost, completePath, http.StatusServiceUnavailable, map[string]string{"code": "nope"})
	putter := &fakePutter{failOn: func(url string) error {
		if strings.HasSuffix(url, "/save/0") {
			return &client.APIError{Method: http.MethodPut, StatusCode: http.StatusForbidden}
		}
		return nil
	}}
	sender := newSaveSender(t, cloud, defaultSaveFiles(t, putter))

	dir := t.TempDir()
	records := []*recordv1.Record{
		saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", writeSaveFile(t, dir, "a.txt", []byte("a")), nil),
		saveRecord(2, savev1.SaveType_SAVE_TYPE_CUSTOM, "b.txt", writeSaveFile(t, dir, "b.txt", []byte("b")), nil),
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 0 || len(result.Remaining) != 2 || len(result.Failures) != 2 {
		t.Fatalf("result = confirmed %d remaining %d failures %d",
			len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
}

// TestSaveCompleteFailure 覆盖 complete 失败：5xx 保留重试，4xx 终结。
func TestSaveCompleteFailure(t *testing.T) {
	cases := []struct {
		name      string
		status    int
		remaining int
		confirmed int
	}{
		{"unavailable", http.StatusServiceUnavailable, 1, 0},
		{"bad request", http.StatusBadRequest, 0, 1},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			serveSavePrepare(cloud)
			cloud.ServeJSON(http.MethodPost, completePath, tc.status, map[string]string{"code": "nope"})
			sender := newSaveSender(t, cloud, defaultSaveFiles(t, &fakePutter{}))

			source := writeSaveFile(t, t.TempDir(), "a.txt", []byte("x"))
			record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", source, nil)
			result := sender.Upload(context.Background(), []*recordv1.Record{record})
			if len(result.Confirmed) != tc.confirmed || len(result.Remaining) != tc.remaining || len(result.Failures) != 1 {
				t.Fatalf("result = confirmed %d remaining %d failures %d",
					len(result.Confirmed), len(result.Remaining), len(result.Failures))
			}
		})
	}
}

// TestSaveInvalidBatchTerminates 覆盖非法分批配置：整批终结且不 prepare。
func TestSaveInvalidBatchTerminates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	files := defaultSaveFiles(t, &fakePutter{})
	files.SaveBatch = 0
	sender := newSaveSender(t, cloud, files)

	source := writeSaveFile(t, t.TempDir(), "a.txt", []byte("x"))
	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", source, nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	if result.Failures[0].Message != saveInvalidBatchPrefix+"0" {
		t.Fatalf("failure = %q", result.Failures[0].Message)
	}
	if got := cloud.CountFor(http.MethodPost, preparePath); got != 0 {
		t.Fatalf("prepare requests = %d, want 0", got)
	}
}

// TestSaveSizeLimitTerminates 覆盖单文件上限：超限直接终结且不 prepare。
func TestSaveSizeLimitTerminates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	files := defaultSaveFiles(t, &fakePutter{})
	files.SaveSize = 3
	sender := newSaveSender(t, cloud, files)

	source := writeSaveFile(t, t.TempDir(), "a.txt", []byte("1234"))
	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "a.txt", source, nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}
	if !strings.HasPrefix(result.Failures[0].Message, saveFileTooLargePrefix) {
		t.Fatalf("failure = %q, want %q prefix", result.Failures[0].Message, saveFileTooLargePrefix)
	}
	if got := cloud.CountFor(http.MethodPost, preparePath); got != 0 {
		t.Fatalf("prepare requests = %d, want 0", got)
	}
}

// TestSaveCustomFallbackUsesWindowsParts 覆盖 CUSTOM 文件按 Windows 语义回退到 files 目录。
func TestSaveCustomFallbackUsesWindowsParts(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveSavePrepare(cloud)
	serveSaveComplete(cloud)
	files := defaultSaveFiles(t, &fakePutter{})
	writeSaveFile(t, files.FilesDir, filepath.Join("sub", "model.py"), []byte("fallback"))
	sender := newSaveSender(t, cloud, files)

	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, `sub\model.py`, "/missing/model.py", nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}
	batches := prepareBatches(t, cloud)
	if len(batches) != 1 || batches[0][0]["path"] != `sub\model.py` {
		t.Fatalf("prepare = %+v, want the original save name", batches)
	}
}

const (
	multipartPreparePath  = "/api/experiment/exp-1/files/prepare-multipart"
	multipartCompletePath = "/api/experiment/exp-1/files/complete-multipart"
)

// serveMultipartPrepare 让假后端按请求里的 count 返回分片 URL。
func serveMultipartPrepare(cloud *clienttest.Cloud) {
	cloud.Handle(http.MethodPost, multipartPreparePath, func(w http.ResponseWriter, _ *http.Request) {
		requests := cloud.Requests()
		var body struct {
			Files []map[string]any `json:"files"`
		}
		_ = json.Unmarshal(requests[len(requests)-1].Body, &body)
		files := make([]map[string]any, 0, len(body.Files))
		for i, item := range body.Files {
			count, _ := item["count"].(float64)
			parts := make([]map[string]any, 0, int(count))
			for n := 1; n <= int(count); n++ {
				parts = append(parts, map[string]any{
					"partNumber": n,
					"url":        fmt.Sprintf("%s/sink/multi/%d/part%d", cloud.URL(), i, n),
				})
			}
			files = append(files, map[string]any{"uploadId": fmt.Sprintf("upload-%d", i), "parts": parts})
		}
		clienttest.WriteJSON(w, http.StatusOK, map[string]any{"files": files})
	})
}

// serveMultipartComplete 让分片合并返回空体（视为全部成功）。
func serveMultipartComplete(cloud *clienttest.Cloud) {
	cloud.ServeJSON(http.MethodPost, multipartCompletePath, http.StatusCreated, nil)
}

// multipartCounts 返回 prepare-multipart 请求里每个文件的分片数。
func multipartCounts(t *testing.T, cloud *clienttest.Cloud) []int {
	t.Helper()
	requests := cloud.RequestsFor(http.MethodPost, multipartPreparePath)
	counts := make([]int, 0)
	for _, request := range requests {
		files, _ := decodeBody(t, request.Body)["files"].([]any)
		for _, item := range files {
			count, _ := item.(map[string]any)["count"].(float64)
			counts = append(counts, int(count))
		}
	}
	return counts
}

// TestSaveSplitRoutesFiles 覆盖 SaveSplit 分界：小于分界走直传，等于分界走分片。
func TestSaveSplitRoutesFiles(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveSavePrepare(cloud)
	serveSaveComplete(cloud)
	serveMultipartPrepare(cloud)
	serveMultipartComplete(cloud)
	putter := &fakePutter{}
	files := defaultSaveFiles(t, putter)
	files.SaveSplit = 10
	files.SavePart = 4
	sender := newSaveSender(t, cloud, files)

	dir := t.TempDir()
	records := []*recordv1.Record{
		saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "small.txt", writeSaveFile(t, dir, "small.txt", []byte("123456789")), nil),
		saveRecord(2, savev1.SaveType_SAVE_TYPE_CUSTOM, "big.bin", writeSaveFile(t, dir, "big.bin", []byte("0123456789")), nil),
	}
	result := sender.Upload(context.Background(), []*recordv1.Record{records[0], records[1]})
	if len(result.Confirmed) != 2 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}

	batches := prepareBatches(t, cloud)
	if len(batches) != 1 || len(batches[0]) != 1 || batches[0][0]["path"] != "small.txt" {
		t.Fatalf("small prepare = %+v, want only small.txt", batches)
	}
	if counts := multipartCounts(t, cloud); len(counts) != 1 || counts[0] != 3 {
		t.Fatalf("multipart counts = %v, want [3]", counts)
	}
	puts := putter.snapshot()
	if len(puts) != 1+3 {
		t.Fatalf("puts = %d, want 1 small + 3 parts", len(puts))
	}
}

// TestSaveMultipartPartBoundaries 覆盖分片数计算与分片内容拼接。
func TestSaveMultipartPartBoundaries(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveMultipartPrepare(cloud)
	serveMultipartComplete(cloud)
	putter := &fakePutter{}
	files := defaultSaveFiles(t, putter)
	files.SaveSplit = 1
	files.SavePart = 4
	sender := newSaveSender(t, cloud, files)

	contents := [][]byte{[]byte("abcdefgh"), []byte("abcdefghi"), []byte("abc")}
	dir := t.TempDir()
	records := make([]*recordv1.Record, 0, len(contents))
	for i, content := range contents {
		name := fmt.Sprintf("f%d.bin", i)
		records = append(records, saveRecord(int64(i+1), savev1.SaveType_SAVE_TYPE_CUSTOM, name, writeSaveFile(t, dir, name, content), nil))
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 3 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}
	if counts := multipartCounts(t, cloud); len(counts) != 3 || counts[0] != 2 || counts[1] != 3 || counts[2] != 1 {
		t.Fatalf("multipart counts = %v, want [2 3 1]", counts)
	}

	byURL := make(map[string][]byte)
	for _, put := range putter.snapshot() {
		byURL[put.url] = put.body
	}
	for i, content := range contents {
		var joined []byte
		count := (len(content) + 3) / 4
		for n := 1; n <= count; n++ {
			url := fmt.Sprintf("%s/sink/multi/%d/part%d", cloud.URL(), i, n)
			part, ok := byURL[url]
			if !ok {
				t.Fatalf("missing part %s", url)
			}
			joined = append(joined, part...)
		}
		if string(joined) != string(content) {
			t.Fatalf("file %d reassembled = %q, want %q", i, joined, content)
		}
	}

	completed := cloud.RequestsFor(http.MethodPost, multipartCompletePath)
	if len(completed) != 1 {
		t.Fatalf("complete-multipart requests = %d, want 1", len(completed))
	}
	filesBody, _ := decodeBody(t, completed[0].Body)["files"].([]any)
	if len(filesBody) != 3 {
		t.Fatalf("complete files = %d, want 3", len(filesBody))
	}
	parts, _ := filesBody[1].(map[string]any)["parts"].([]any)
	if len(parts) != 3 {
		t.Fatalf("second file parts = %d, want 3", len(parts))
	}
}

// readThenFailPutter 读完整片后按 URL 返回失败，用于观测进度回滚。
type readThenFailPutter struct {
	mu      sync.Mutex
	puts    []objectPut
	failURL string
	err     error
}

// PutObject 实现 ObjectPutter。
func (p *readThenFailPutter) PutObject(_ context.Context, url string, body io.Reader, _ int64, contentType string) (string, error) {
	data, err := io.ReadAll(body)
	if err != nil {
		return "", err
	}
	p.mu.Lock()
	p.puts = append(p.puts, objectPut{url: url, contentType: contentType, body: data})
	p.mu.Unlock()
	if strings.HasSuffix(url, p.failURL) {
		return "", p.err
	}
	return "etag", nil
}

// TestSaveMultipartPartRetryableRollsBackProgress 覆盖单片失败：记录保留重试，
// 该分片已上报的进度回滚到 0。
func TestSaveMultipartPartRetryableRollsBackProgress(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveMultipartPrepare(cloud)
	putter := &readThenFailPutter{
		failURL: "/part2",
		err:     &client.APIError{Method: http.MethodPut, StatusCode: http.StatusInternalServerError},
	}
	tracker := uploader.NewTracker()
	files := defaultSaveFiles(t, putter)
	files.SaveSplit = 1
	files.SavePart = 4
	files.Tracker = tracker
	sender := newSaveSender(t, cloud, files)

	dir := t.TempDir()
	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "big.bin", writeSaveFile(t, dir, "big.bin", []byte("0123456789")), nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 0 || len(result.Remaining) != 1 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	if got := cloud.CountFor(http.MethodPost, multipartCompletePath); got != 0 {
		t.Fatalf("complete-multipart requests = %d, want 0", got)
	}
	stats := tracker.Snapshot()
	if stats.TotalNumber != 1 || stats.UploadedNumber != 0 || len(stats.Files) != 1 {
		t.Fatalf("tracker = total %d uploaded %d files %d", stats.TotalNumber, stats.UploadedNumber, len(stats.Files))
	}
	// 分片 1、3 成功（4 + 2 字节），分片 2 读完 4 字节后失败并回滚，总进度应为 6。
	if got := stats.Files[0].Uploaded; got != 6 {
		t.Fatalf("uploaded = %d, want 6 after rolling back the failed part", got)
	}
}

// TestSaveMultipartLocalErrorTerminates 覆盖分片阶段本地读错：终结并上报 FAILED。
func TestSaveMultipartLocalErrorTerminates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	serveMultipartPrepare(cloud)
	serveSaveComplete(cloud)
	putter := &fakePutter{failOn: func(url string) error {
		if strings.HasSuffix(url, "/part1") {
			return &fs.PathError{Op: "open", Path: "big.bin", Err: fs.ErrNotExist}
		}
		return nil
	}}
	tracker := uploader.NewTracker()
	files := defaultSaveFiles(t, putter)
	files.SaveSplit = 1
	files.SavePart = 4
	files.Tracker = tracker
	sender := newSaveSender(t, cloud, files)

	dir := t.TempDir()
	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "big.bin", writeSaveFile(t, dir, "big.bin", []byte("0123456789")), nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	completed := completeBatches(t, cloud)
	if len(completed) != 1 || len(completed[0]) != 1 || completed[0][0]["state"] != saveStateFailed {
		t.Fatalf("complete = %+v, want one FAILED entry", completed)
	}
	if got := cloud.CountFor(http.MethodPost, multipartCompletePath); got != 0 {
		t.Fatalf("complete-multipart requests = %d, want 0", got)
	}
	if stats := tracker.Snapshot(); stats.TotalNumber != 1 || stats.UploadedNumber != 0 || len(stats.Files) != 0 {
		t.Fatalf("tracker = total %d uploaded %d files %d", stats.TotalNumber, stats.UploadedNumber, len(stats.Files))
	}
}

// TestSaveMultipartCompleteStates 覆盖 complete-multipart 的三种响应。
func TestSaveMultipartCompleteStates(t *testing.T) {
	cases := []struct {
		name      string
		status    int
		payload   map[string]any
		confirmed int
		remaining int
		failures  int
	}{
		{"uploaded", http.StatusCreated, nil, 1, 0, 0},
		{"failed state", http.StatusCreated, map[string]any{
			"files": []any{map[string]any{"path": "big.bin", "state": "FAILED"}},
		}, 1, 0, 1},
		{"unavailable", http.StatusServiceUnavailable, map[string]any{"code": "nope"}, 0, 1, 1},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			serveMultipartPrepare(cloud)
			if tc.status == http.StatusServiceUnavailable {
				cloud.ServeJSON(http.MethodPost, multipartCompletePath, tc.status, tc.payload)
			} else if tc.payload == nil {
				serveMultipartComplete(cloud)
			} else {
				cloud.ServeJSON(http.MethodPost, multipartCompletePath, tc.status, tc.payload)
			}
			files := defaultSaveFiles(t, &fakePutter{})
			files.SaveSplit = 1
			files.SavePart = 4
			sender := newSaveSender(t, cloud, files)

			dir := t.TempDir()
			record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "big.bin", writeSaveFile(t, dir, "big.bin", []byte("0123456789")), nil)
			result := sender.Upload(context.Background(), []*recordv1.Record{record})
			if len(result.Confirmed) != tc.confirmed || len(result.Remaining) != tc.remaining || len(result.Failures) != tc.failures {
				t.Fatalf("result = confirmed %d remaining %d failures %d",
					len(result.Confirmed), len(result.Remaining), len(result.Failures))
			}
			if tc.failures > 0 && !strings.HasPrefix(result.Failures[0].Message, saveMultipartStateFailedPrefix) &&
				tc.status == http.StatusCreated {
				t.Fatalf("failure = %q, want %q prefix", result.Failures[0].Message, saveMultipartStateFailedPrefix)
			}
		})
	}
}

// TestSaveMultipartETagForwarded 覆盖分片 ETag 原样上报且直传不带鉴权头。
func TestSaveMultipartETagForwarded(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodGet, "/api/auth/verify", http.StatusOK, map[string]any{
		"uid": 42, "username": "alice", "createdAt": "2026-01-01T00:00:00.000Z",
	})
	serveMultipartPrepare(cloud)
	serveMultipartComplete(cloud)
	cloud.Handle(http.MethodPut, "/sink/multi/0/part1", func(w http.ResponseWriter, _ *http.Request) {
		w.Header().Set("ETag", `"quoted-etag"`)
		w.WriteHeader(http.StatusOK)
	})

	transport, err := client.New(context.Background(), cloud.URL(), "smoke-key", client.Options{})
	if err != nil {
		t.Fatalf("create client: %v", err)
	}
	files := defaultSaveFiles(t, transport)
	files.SaveSplit = 1
	files.SavePart = 10
	sender := NewSave(Options{
		API:          api.New(transport),
		Username:     "alice",
		Project:      "demo",
		ProjectID:    "proj-1",
		ExperimentID: "exp-1",
		Files:        files,
	})

	dir := t.TempDir()
	record := saveRecord(1, savev1.SaveType_SAVE_TYPE_CUSTOM, "big.bin", writeSaveFile(t, dir, "big.bin", []byte("12345")), nil)
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}

	puts := cloud.RequestsFor(http.MethodPut, "/sink/multi/0/part1")
	if len(puts) != 1 {
		t.Fatalf("part requests = %d, want 1", len(puts))
	}
	if got := puts[0].Header.Get("Authorization"); got != "" {
		t.Fatalf("part Authorization = %q, want none", got)
	}
	completed := cloud.RequestsFor(http.MethodPost, multipartCompletePath)
	if len(completed) != 1 {
		t.Fatalf("complete-multipart requests = %d, want 1", len(completed))
	}
	filesBody, _ := decodeBody(t, completed[0].Body)["files"].([]any)
	parts, _ := filesBody[0].(map[string]any)["parts"].([]any)
	if len(parts) != 1 || parts[0].(map[string]any)["etag"] != "quoted-etag" {
		t.Fatalf("parts = %+v, want the ETag without quotes", parts)
	}
}
