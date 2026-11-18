// Package runctx 维护单个 run 的运行上下文：配置、本地目录布局与云端在线参数。
//
// 配置字段由内嵌的 config.Settings 提供，New 之后视为只读。路径方法只做
// 字符串拼接，不创建目录。打印配置时必须走 config.MaskKey，避免泄漏 API Key。
package runctx

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"slices"

	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/fs"
)

// runIDMaxLen 是 run 文件名中 run_id 的最大字节长度。
const runIDMaxLen = 100

var (
	// ErrNoSyncRunFile 表示 run_dir 中不存在 run-*.swanlab 文件。
	ErrNoSyncRunFile = errors.New("no run-*.swanlab file found")
	// ErrMultipleSyncRunFiles 表示 run_dir 中存在多个 run-*.swanlab 文件。
	ErrMultipleSyncRunFiles = errors.New("multiple run-*.swanlab files found")
)

// RunContext 是单个 run 的运行上下文。
type RunContext struct {
	config.Settings

	online    OnlineParams
	hasOnline bool
}

// OnlineParams 是云端身份信息，由 prepare_experiment_start 的响应注入。
type OnlineParams struct {
	Username       string
	Project        string
	ProjectID      string
	ProjectVersion int32
	ExperimentID   string
}

// New 校验配置并创建运行上下文。
func New(s config.Settings) (*RunContext, error) {
	if err := s.Validate(); err != nil {
		return nil, err
	}
	return &RunContext{Settings: s}, nil
}

// MediaDir 返回媒体文件目录：run_dir/media。
func (c *RunContext) MediaDir() string { return filepath.Join(c.RunDir, "media") }

// DebugDir 返回诊断日志目录：run_dir/debug。
func (c *RunContext) DebugDir() string { return filepath.Join(c.RunDir, "debug") }

// FilesDir 返回保存文件目录：run_dir/files。
func (c *RunContext) FilesDir() string { return filepath.Join(c.RunDir, "files") }

// MetadataFile 返回内部元数据文件路径：files/swanlab-metadata.json。
func (c *RunContext) MetadataFile() string {
	return filepath.Join(c.FilesDir(), "swanlab-metadata.json")
}

// ConfigFile 返回用户配置文件路径：files/config.yaml。
func (c *RunContext) ConfigFile() string { return filepath.Join(c.FilesDir(), "config.yaml") }

// RequirementsFile 返回依赖文件路径：files/requirements.txt。
func (c *RunContext) RequirementsFile() string {
	return filepath.Join(c.FilesDir(), "requirements.txt")
}

// CondaFile 返回 conda 环境文件路径：files/conda.yaml。
func (c *RunContext) CondaFile() string { return filepath.Join(c.FilesDir(), "conda.yaml") }

// RunFile 返回本地日志文件路径：run_dir/run-{run_id≤100}.swanlab。
// run_id 先经 fs.SafeFmt 清洗，再截断到 100 字节。
func (c *RunContext) RunFile() (string, error) {
	if c.RunID == "" {
		return "", fmt.Errorf("runctx: %w", config.ErrRunIDRequired)
	}
	id, _, err := fs.SafeTruncate(fs.SafeFmt(c.RunID, "unknown_run_id"), runIDMaxLen)
	if err != nil {
		return "", fmt.Errorf("runctx: derive run file for run_id %q: %w", c.RunID, err)
	}
	return filepath.Join(c.RunDir, "run-"+id+".swanlab"), nil
}

// CheckpointFile 返回上传断点文件路径：RunFile + ".checkpoint.json"。
func (c *RunContext) CheckpointFile() (string, error) {
	runFile, err := c.RunFile()
	if err != nil {
		return "", err
	}
	return runFile + ".checkpoint.json", nil
}

// SetOnlineParams 注入云端身份信息，重复调用以最后一次为准。
func (c *RunContext) SetOnlineParams(p OnlineParams) {
	c.online = p
	c.hasOnline = true
}

// Online 返回已注入的云端身份信息；未注入时第二个返回值为 false。
func (c *RunContext) Online() (OnlineParams, bool) {
	return c.online, c.hasOnline
}

// FindSyncRunFile 在 runDir 下查找唯一的 run-*.swanlab 文件供 sync 读取。
// sync 模式下目标 run_id 可能被覆盖或新生成，因此不能按 RunContext.RunID 拼接文件名。
func FindSyncRunFile(runDir string) (string, error) {
	if runDir == "" {
		return "", fmt.Errorf("runctx: %w", config.ErrRunDirRequired)
	}
	info, err := os.Stat(runDir)
	if err != nil {
		return "", fmt.Errorf("runctx: stat run_dir %q: %w", runDir, err)
	}
	if !info.IsDir() {
		return "", fmt.Errorf("runctx: run_dir %q is not a directory", runDir)
	}
	matches, err := filepath.Glob(filepath.Join(runDir, "run-*.swanlab"))
	if err != nil {
		return "", fmt.Errorf("runctx: scan run files in %q: %w", runDir, err)
	}
	slices.Sort(matches)
	switch len(matches) {
	case 0:
		return "", fmt.Errorf("runctx: %w in %s", ErrNoSyncRunFile, runDir)
	case 1:
		fileInfo, err := os.Stat(matches[0])
		if err != nil {
			return "", fmt.Errorf("runctx: stat run file %q: %w", matches[0], err)
		}
		if fileInfo.IsDir() {
			return "", fmt.Errorf("runctx: run file %q is a directory", matches[0])
		}
		return matches[0], nil
	default:
		return "", fmt.Errorf("runctx: %w in %s", ErrMultipleSyncRunFiles, runDir)
	}
}
