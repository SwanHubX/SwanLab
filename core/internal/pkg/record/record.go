// Package record 冻结数据面内部共享的记录口径。
//
// 记录本体直接复用 proto 的 swanlab.record.v1.Record，不另建平行结构：
// 同一份消息既用于生成落盘字节，也用于编码上传请求，避免两套字段各自演化。
//
// 本包只固定所有下游必须一致的三件事：
//   - Type：oneof record_type 的稳定标签，intake 分流与 uploader 分轨都用它；
//   - TypeOf：从消息推导标签，保证标签与 oneof 不会漂移；
//   - Size：配额计量口径，intake 预留与 uploader 归还必须调用同一个函数。
package record

import (
	"google.golang.org/protobuf/proto"

	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// Type 是 Record 的 oneof 分支标签，取值与 proto 字段名一致。
type Type string

const (
	TypeStart  Type = "start"
	TypeFinish Type = "finish"
	TypeColumn Type = "column"
	TypeScalar Type = "scalar"
	TypeMedia  Type = "media"
	TypeLog    Type = "log"
	TypeSave   Type = "save"
)

// 固定槽位号，见 record.proto 对 Record.num 的约定：
// Start 与 Finish 仅写本地，其余为可重复上报的 Save 记录。
// 正数编号从 1 开始由 builder 递增分配，小于 -100 的区间留给离线同步。
const (
	NumStart            int64 = -1
	NumFinish           int64 = -2
	NumConfigSave       int64 = -3
	NumMetadataSave     int64 = -4
	NumRequirementsSave int64 = -5
	NumCondaSave        int64 = -6
)

// TypeOf 返回 r 的 oneof 标签；r 为 nil 或 oneof 未设置时返回空字符串。
func TypeOf(r *recordv1.Record) Type {
	switch r.GetRecordType().(type) {
	case *recordv1.Record_Start:
		return TypeStart
	case *recordv1.Record_Finish:
		return TypeFinish
	case *recordv1.Record_Column:
		return TypeColumn
	case *recordv1.Record_Scalar:
		return TypeScalar
	case *recordv1.Record_Media:
		return TypeMedia
	case *recordv1.Record_Log:
		return TypeLog
	case *recordv1.Record_Save:
		return TypeSave
	default:
		return ""
	}
}

// Size 返回记录序列化后的字节数（不含 store 的 7 字节物理块头），
// 是 store 通道容量与 UploadBudget 的唯一计量口径。
func Size(r *recordv1.Record) int64 {
	if r == nil {
		return 0
	}
	return int64(proto.Size(r))
}
