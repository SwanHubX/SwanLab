// Package runctx 维护单个 run 的运行上下文：配置、本地目录布局与云端在线参数。
//
// 配置字段由内嵌的 config.Settings 提供，New 之后视为只读。路径方法只做
// 字符串拼接，不创建目录。打印配置时必须走 config.MaskKey，避免泄漏 API Key。
package runctx

import (
	"fmt"
	"path/filepath"

	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/pkg/fs"
)

// runIDMaxLen 是 run 文件名中 run_id 的最大字节长度。
const runIDMaxLen = 100

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

// SetOnlineParams 注入云端身份信息，重复调用以最后一次为准。
func (c *RunContext) SetOnlineParams(p OnlineParams) {
	c.online = p
	c.hasOnline = true
}

// OnlineParams 返回已注入的云端身份信息；未注入时第二个返回值为 false。
func (c *RunContext) OnlineParams() (OnlineParams, bool) {
	return c.online, c.hasOnline
}
