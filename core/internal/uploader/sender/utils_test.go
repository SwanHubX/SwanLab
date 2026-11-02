package sender

import (
	"context"
	"errors"
	"io/fs"
	"net/http"
	"path/filepath"
	"testing"

	"golang.org/x/sync/semaphore"

	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

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

// TestPutErrorRetryable 覆盖对象存储 PUT 失败的分类。
func TestPutErrorRetryable(t *testing.T) {
	cases := []struct {
		name string
		err  error
		want bool
	}{
		{"forbidden", &client.APIError{Method: http.MethodPut, StatusCode: http.StatusForbidden}, true},
		{"server error", &client.APIError{Method: http.MethodPut, StatusCode: http.StatusInternalServerError}, true},
		{"transport", &client.APIError{Method: http.MethodPut, Code: "network_error", Retryable: true}, true},
		{"protocol", client.NewProtocolError(http.MethodPut, "https://s3", "bad_response", "bad", nil), false},
		{"local", &fs.PathError{Op: "open", Path: "a.txt", Err: fs.ErrNotExist}, false},
		{"canceled", context.Canceled, true},
		{"deadline", context.DeadlineExceeded, true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := putErrorRetryable(tc.err); got != tc.want {
				t.Fatalf("putErrorRetryable(%v) = %v, want %v", tc.err, got, tc.want)
			}
		})
	}
}

// TestComputeMD5s 覆盖批量 MD5 的顺序与错误槽位。
func TestComputeMD5s(t *testing.T) {
	dir := t.TempDir()
	paths := []string{
		writeSaveFile(t, dir, "a.txt", []byte("abc")),
		filepath.Join(dir, "missing"),
	}
	digests, errs := computeMD5s(context.Background(), paths)
	if digests[0] != "900150983cd24fb0d6963f7d28e17f72" {
		t.Fatalf("digest = %q", digests[0])
	}
	if errs[0] != nil {
		t.Fatalf("errs[0] = %v", errs[0])
	}
	if errs[1] == nil {
		t.Fatal("missing file should produce an error")
	}
}

// TestMergeResult 覆盖三态结果的合并。
func TestMergeResult(t *testing.T) {
	first := recordv1.Record{Num: 1}
	second := recordv1.Record{Num: 2}
	merged := mergeResult(
		uploader.Result{Confirmed: []*recordv1.Record{&first}, Failures: []uploader.Failure{{Message: "a"}}},
		uploader.Result{Remaining: []*recordv1.Record{&second}, Failures: []uploader.Failure{{Message: "b"}}},
	)
	if len(merged.Confirmed) != 1 || merged.Confirmed[0] != &first {
		t.Fatalf("confirmed = %v", merged.Confirmed)
	}
	if len(merged.Remaining) != 1 || merged.Remaining[0] != &second {
		t.Fatalf("remaining = %v", merged.Remaining)
	}
	if len(merged.Failures) != 2 {
		t.Fatalf("failures = %d, want 2", len(merged.Failures))
	}
}

// TestRunConcurrently 覆盖错误按索引落位与 ctx 取消。
func TestRunConcurrently(t *testing.T) {
	errs := runConcurrently(context.Background(), semaphore.NewWeighted(2), 4, func(i int) error {
		if i%2 == 1 {
			return errors.New("odd")
		}
		return nil
	})
	for i, err := range errs {
		if (i%2 == 1) != (err != nil) {
			t.Fatalf("errs[%d] = %v", i, err)
		}
	}

	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	canceled := runConcurrently(ctx, semaphore.NewWeighted(2), 3, func(int) error { return nil })
	for i, err := range canceled {
		if !errors.Is(err, context.Canceled) {
			t.Fatalf("errs[%d] = %v, want context.Canceled", i, err)
		}
	}
}
