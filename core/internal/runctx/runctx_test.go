package runctx

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/config"
)

// newContext 用默认配置构造运行上下文，只覆盖测试关心的字段。
func newContext(t *testing.T, runID, runDir string) *RunContext {
	t.Helper()
	cfg := config.Defaults()
	cfg.RunID, cfg.RunDir = runID, runDir
	c, err := New(cfg)
	if err != nil {
		t.Fatalf("New(%+v) 报错: %v", cfg, err)
	}
	return c
}

// newEmptyContext 绕过 New 的校验，构造缺少 run_id 的上下文，
// 用于验证 RunFile 自身的兜底错误。
func newEmptyContext() *RunContext {
	cfg := config.Defaults()
	cfg.RunDir = "/rd"
	return &RunContext{Settings: cfg}
}

func TestContextPaths(t *testing.T) {
	c := newContext(t, "r1", filepath.Join("tmp", "run"))
	if got := c.MediaDir(); got != filepath.Join("tmp", "run", "media") {
		t.Fatalf("MediaDir = %q", got)
	}
	if got := c.FilesDir(); got != filepath.Join("tmp", "run", "files") {
		t.Fatalf("FilesDir = %q", got)
	}
}

func TestRunFile(t *testing.T) {
	tests := []struct {
		name, runID, wantBase string
	}{
		{"常规", "abc123", "run-abc123.swanlab"},
		{"特殊字符清洗", "a<b>c:d\"e/f", "run-a_b_c_d_e_f.swanlab"},
		{"点横线替换", "exp.1-x", "run-exp_1_x.swanlab"},
		{"超长截断", strings.Repeat("id", 80), "run-" + "idi" + "..." + strings.Repeat("id", 47) + ".swanlab"},
		{"全非法回退", "???###", "run-unknown_run_id.swanlab"},
	}
	runDir := filepath.Join("tmp", "run")
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			c := newContext(t, tc.runID, runDir)
			got, err := c.RunFile()
			if err != nil {
				t.Fatalf("RunFile() 报错: %v", err)
			}
			if filepath.Base(got) != tc.wantBase {
				t.Fatalf("RunFile() base = %q, want %q", filepath.Base(got), tc.wantBase)
			}
			if filepath.Dir(got) != runDir {
				t.Fatalf("RunFile() dir = %q, want %q", filepath.Dir(got), runDir)
			}
			// 文件名字节长度约束：run-(≤100).swanlab
			if n := len(filepath.Base(got)); n > len("run-.swanlab")+100 {
				t.Fatalf("文件名 %d 字节超过上限", n)
			}
		})
	}

	if _, err := newEmptyContext().RunFile(); !errors.Is(err, config.ErrRunIDRequired) {
		t.Fatalf("空 run_id 应报错, got %v", err)
	}
}

func TestNewValidates(t *testing.T) {
	cfg := config.Defaults()
	cfg.RunDir = "/rd"
	if _, err := New(cfg); !errors.Is(err, config.ErrRunIDRequired) {
		t.Fatalf("缺 run_id: %v", err)
	}

	cfg = config.Defaults()
	cfg.RunID = "r"
	if _, err := New(cfg); !errors.Is(err, config.ErrRunDirRequired) {
		t.Fatalf("缺 run_dir: %v", err)
	}

	// 数值边界校验来自 config，runctx 原样透传。
	cfg = config.Defaults()
	cfg.RunID, cfg.RunDir = "r", "/rd"
	cfg.RecordBatch = 0
	if _, err := New(cfg); err == nil {
		t.Fatal("record_batch=0 应被拒绝")
	}

	// 合法配置原样保留字段。
	cfg = config.Defaults()
	cfg.RunID, cfg.RunDir, cfg.SectionRule = "r", "/rd", 1
	cfg.APIKey, cfg.APIHost = "0123456789abcdef", "https://api.example.com"
	c, err := New(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if c.RunID != "r" || c.SectionRule != 1 || c.APIKey != cfg.APIKey || c.APIHost != cfg.APIHost {
		t.Fatalf("配置字段丢失: %+v", c)
	}
	// 内嵌的 String() 被提升，且默认走脱敏口径。
	if out := c.String(); strings.Contains(out, cfg.APIKey) || !strings.Contains(out, config.MaskKey(cfg.APIKey)) {
		t.Fatalf("String() 未脱敏: %s", out)
	}
}

func TestOnlineParams(t *testing.T) {
	c := newContext(t, "r", "/rd")
	if _, ok := c.Online(); ok {
		t.Fatal("未注入时 Online 应返回 false")
	}
	p := OnlineParams{
		Username:     "u",
		Project:      "p",
		ProjectID:    "pid",
		ExperimentID: "eid",
	}
	c.SetOnlineParams(p)
	got, ok := c.Online()
	if !ok || got != p {
		t.Fatalf("Online = (%+v, %v)", got, ok)
	}
}

func TestOnlineParamsConcurrent(t *testing.T) {
	c := newContext(t, "r", "/rd")
	const writers = 4
	const readers = 8
	const iters = 200

	var wg sync.WaitGroup
	for w := range writers {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			for i := range iters {
				tag := string(rune('a' + (id+i)%26))
				c.SetOnlineParams(OnlineParams{
					Username:       "u-" + tag,
					Project:        "p-" + tag,
					ProjectID:      "pid-" + tag,
					ProjectVersion: int32(id + 1),
					ExperimentID:   "eid-" + tag,
				})
			}
		}(w)
	}
	for range readers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for range iters {
				got, ok := c.Online()
				if ok && (got.Username == "" || got.ProjectID == "" || got.ExperimentID == "") {
					t.Errorf("observed partial OnlineParams: %+v", got)
					return
				}
			}
		}()
	}
	wg.Wait()

	got, ok := c.Online()
	if !ok || got.Username == "" || got.ProjectID == "" || got.ExperimentID == "" {
		t.Fatalf("final Online = (%+v, %v)", got, ok)
	}
}

func TestCheckpointFile(t *testing.T) {
	c := newContext(t, "run-123", "/rd")
	cp, err := c.CheckpointFile()
	if err != nil {
		t.Fatal(err)
	}
	runFile, err := c.RunFile()
	if err != nil {
		t.Fatal(err)
	}
	if cp != runFile+".checkpoint.json" {
		t.Fatalf("CheckpointFile = %q, want %q", cp, runFile+".checkpoint.json")
	}

	c.RunID = ""
	if _, err := c.CheckpointFile(); !errors.Is(err, config.ErrRunIDRequired) {
		t.Fatalf("空 RunID 应返回 ErrRunIDRequired, got %v", err)
	}
}

func TestFindSyncRunFile(t *testing.T) {
	t.Run("single file", func(t *testing.T) {
		dir := t.TempDir()
		want := filepath.Join(dir, "run-abc.swanlab")
		if err := os.WriteFile(want, nil, 0o644); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(dir, "other.txt"), nil, 0o644); err != nil {
			t.Fatal(err)
		}
		got, err := FindSyncRunFile(dir)
		if err != nil {
			t.Fatalf("FindSyncRunFile: %v", err)
		}
		if got != want {
			t.Fatalf("FindSyncRunFile = %q, want %q", got, want)
		}
	})

	t.Run("missing file", func(t *testing.T) {
		dir := t.TempDir()
		if _, err := FindSyncRunFile(dir); !errors.Is(err, ErrNoSyncRunFile) {
			t.Fatalf("got %v, want ErrNoSyncRunFile", err)
		}
	})

	t.Run("multiple files", func(t *testing.T) {
		dir := t.TempDir()
		for _, name := range []string{"run-a.swanlab", "run-b.swanlab"} {
			if err := os.WriteFile(filepath.Join(dir, name), nil, 0o644); err != nil {
				t.Fatal(err)
			}
		}
		if _, err := FindSyncRunFile(dir); !errors.Is(err, ErrMultipleSyncRunFiles) {
			t.Fatalf("got %v, want ErrMultipleSyncRunFiles", err)
		}
	})

	t.Run("invalid run_dir", func(t *testing.T) {
		if _, err := FindSyncRunFile(""); !errors.Is(err, config.ErrRunDirRequired) {
			t.Fatalf("empty dir got %v, want ErrRunDirRequired", err)
		}
		file := filepath.Join(t.TempDir(), "not-a-dir")
		if err := os.WriteFile(file, nil, 0o644); err != nil {
			t.Fatal(err)
		}
		if _, err := FindSyncRunFile(file); err == nil {
			t.Fatal("non-directory run_dir should fail")
		}
		if _, err := FindSyncRunFile(filepath.Join(t.TempDir(), "missing")); !errors.Is(err, os.ErrNotExist) {
			t.Fatalf("missing dir got %v, want ErrNotExist", err)
		}
	})

	t.Run("directory named like run file", func(t *testing.T) {
		dir := t.TempDir()
		if err := os.Mkdir(filepath.Join(dir, "run-dir.swanlab"), 0o755); err != nil {
			t.Fatal(err)
		}
		if _, err := FindSyncRunFile(dir); err == nil {
			t.Fatal("directory matching run-*.swanlab should fail")
		}
	})
}
