// Package config 定义单个 run 的运行配置。
//
// Settings 与 CoreSettings 协议一一对应，是纯 Go 数据：不 import proto，
// 也不依赖任何 internal 业务包。本包只负责字段、默认值、边界校验与日志
// 脱敏；目录路径、在线参数等运行期状态属于 runctx。
//
// 环境变量与配置文件由客户端解析，Go 只接收解析后的值。
package config

import (
	"errors"
	"fmt"
	"math"
)

// 默认值与 Python 客户端一致，供测试、工具与冒烟驱动使用。
// 客户端总是显式下发这些字段，Validate 不会静默补零值。
const (
	// DefaultRecordBatch 是单批上传的记录条数上限。
	DefaultRecordBatch = 10_000
	// DefaultRecordInterval 是攒批的最大等待秒数。
	DefaultRecordInterval = 1.5
	// DefaultSaveSize 是单文件保存上限（50 GiB）。
	DefaultSaveSize int64 = 50 << 30
	// DefaultSaveSplit 是大小文件分界（100 MiB）。
	DefaultSaveSplit int64 = 100 << 20
	// DefaultSavePart 是大文件分片大小（32 MiB）。
	DefaultSavePart int64 = 32 << 20
	// DefaultSaveBatch 是小文件上传的单批数量。
	DefaultSaveBatch = 100
)

// recordBatchMax 是 record_batch 的上界（开区间）。
const recordBatchMax = 100_000

// maskEdge 是脱敏时保留的首尾字节数；maskMinLen 是启用首尾保留的最短长度。
const (
	maskEdge   = 4
	maskMinLen = maskEdge * 2
)

// 必填项哨兵错误，供调用方与测试用 errors.Is 判别。
var (
	// ErrRunIDRequired 表示缺少 run_id。
	ErrRunIDRequired = errors.New("run_id is required")
	// ErrRunDirRequired 表示缺少 run_dir。
	ErrRunDirRequired = errors.New("run_dir is required")
)

// Settings 是单个会话的运行配置。
//
// 字段在构造后视为只读；云端身份（username/project/experiment_id 等）属于
// 运行期状态，由 runctx 持有，不放在这里。
type Settings struct {
	// RunID 会话唯一标识，同时决定本地 run-{id}.swanlab 文件名。
	RunID string
	// RunDir 本次运行的本地根目录（绝对路径，OS 相关格式）。
	RunDir string

	// SectionRule 自动列的 section 切分规则（见 builder.BuildAutoColumn）。
	SectionRule int32
	// RecordBatch 单批上传的记录条数上限。
	RecordBatch int32
	// RecordInterval 攒批的最大等待秒数。
	RecordInterval float64
	// SaveSize 单文件保存上限（字节）。
	SaveSize int64
	// SaveSplit 大小文件分界（字节）。
	SaveSplit int64
	// SavePart 大文件分片大小（字节）。
	SavePart int64
	// SaveBatch 小文件上传的单批数量。
	SaveBatch int32

	// SkipStore 跳过本地 Record 落盘（仅 online 模式合法，由上游校验）。
	// 启动期配置，运行期不得修改；磁盘写满等故障属于会话运行期状态。
	SkipStore bool

	// APIKey 云端认证凭据；日志输出必须脱敏。
	APIKey string
	// APIHost 云端 API 站点（如 https://api.swanlab.cn）。
	APIHost string
	// WebHost 云端 Web 站点。
	WebHost string
}

// Defaults 返回与 Python CoreSettings 一致的默认配置。
//
// 只填有非零默认值的字段：RunID / RunDir / SectionRule / SkipStore /
// APIKey / APIHost / WebHost 保持零值（SectionRule 的默认值就是 0）。
func Defaults() Settings {
	return Settings{
		RecordBatch:    DefaultRecordBatch,
		RecordInterval: DefaultRecordInterval,
		SaveSize:       DefaultSaveSize,
		SaveSplit:      DefaultSaveSplit,
		SavePart:       DefaultSavePart,
		SaveBatch:      DefaultSaveBatch,
	}
}

// Validate 校验必填项与数值边界。
//
// record_batch 必须落在 (0, 100_000)；record_interval 与 save_part 必须为正，
// record_interval 额外拒绝 NaN（+Inf 允许）；其余 save_* 不限制。
func (s Settings) Validate() error {
	if s.RunID == "" {
		return fmt.Errorf("config: %w", ErrRunIDRequired)
	}
	if s.RunDir == "" {
		return fmt.Errorf("config: %w", ErrRunDirRequired)
	}
	if s.RecordBatch <= 0 || int64(s.RecordBatch) >= recordBatchMax {
		return fmt.Errorf("config: record_batch (%d) must be in (0, %d)", s.RecordBatch, recordBatchMax)
	}
	if math.IsNaN(s.RecordInterval) || s.RecordInterval <= 0 {
		return fmt.Errorf("config: record_interval (%v) must be positive", s.RecordInterval)
	}
	if s.SavePart <= 0 {
		return fmt.Errorf("config: save_part (%d) must be positive", s.SavePart)
	}
	return nil
}

// MaskKey 对 API Key 做日志脱敏：过短的 key 只显示 "****"，
// 其余保留前 maskEdge 与后 maskEdge 字节，中间以 "..." 代替。
func MaskKey(key string) string {
	if key == "" {
		return ""
	}
	if len(key) <= maskMinLen {
		return "****"
	}
	return key[:maskEdge] + "..." + key[len(key)-maskEdge:]
}

// String 返回可安全打印的配置摘要（fmt.Stringer），
// API Key 一律脱敏，避免配置日志泄漏凭据。
func (s Settings) String() string {
	return fmt.Sprintf(
		"run_id=%s run_dir=%s skip_store=%t section_rule=%d record_batch=%d record_interval=%.1f "+
			"save_size=%d save_split=%d save_part=%d save_batch=%d api_host=%s web_host=%s api_key=%s",
		s.RunID, s.RunDir, s.SkipStore, s.SectionRule, s.RecordBatch, s.RecordInterval,
		s.SaveSize, s.SaveSplit, s.SavePart, s.SaveBatch, s.APIHost, s.WebHost, MaskKey(s.APIKey),
	)
}
