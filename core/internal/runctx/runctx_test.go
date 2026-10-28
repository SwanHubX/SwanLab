package runctx

import (
	"errors"
	"path/filepath"
	"strings"
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
	if got := c.DebugDir(); got != filepath.Join("tmp", "run", "debug") {
		t.Fatalf("DebugDir = %q", got)
	}
	if got := c.FilesDir(); got != filepath.Join("tmp", "run", "files") {
		t.Fatalf("FilesDir = %q", got)
	}
	if got := c.MetadataFile(); got != filepath.Join("tmp", "run", "files", "swanlab-metadata.json") {
		t.Fatalf("MetadataFile = %q", got)
	}
	if got := c.ConfigFile(); got != filepath.Join("tmp", "run", "files", "config.yaml") {
		t.Fatalf("ConfigFile = %q", got)
	}
	if got := c.RequirementsFile(); got != filepath.Join("tmp", "run", "files", "requirements.txt") {
		t.Fatalf("RequirementsFile = %q", got)
	}
	if got := c.CondaFile(); got != filepath.Join("tmp", "run", "files", "conda.yaml") {
		t.Fatalf("CondaFile = %q", got)
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
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			c := newContext(t, tc.runID, "/rd")
			got, err := c.RunFile()
			if err != nil {
				t.Fatalf("RunFile() 报错: %v", err)
			}
			if filepath.Base(got) != tc.wantBase {
				t.Fatalf("RunFile() base = %q, want %q", filepath.Base(got), tc.wantBase)
			}
			if filepath.Dir(got) != "/rd" {
				t.Fatalf("RunFile() dir = %q", filepath.Dir(got))
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
	if _, ok := c.OnlineParams(); ok {
		t.Fatal("未注入时 OnlineParams 应返回 false")
	}
	p := OnlineParams{
		Username:     "u",
		Project:      "p",
		ProjectID:    "pid",
		ExperimentID: "eid",
	}
	c.SetOnlineParams(p)
	got, ok := c.OnlineParams()
	if !ok || got != p {
		t.Fatalf("OnlineParams = (%+v, %v)", got, ok)
	}
}
