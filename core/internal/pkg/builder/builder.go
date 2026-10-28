// Package builder 构建数据面 Record 与列定义，只做无状态拼装：
//   - 固定槽位与正数递增编号复用 internal/pkg/record 的常量；
//   - 自动列（Section 切分）、Resume 降级列、默认实验名/颜色词表与 Python 客户端一致；
//   - 所有函数不修改入参消息，时间戳取构建时刻。
//
// 正数编号从 1 开始由 counter 递增分配；小于 -100 的负数区间留给离线同步。
package builder

import (
	"fmt"

	"google.golang.org/protobuf/proto"
	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

// BuildColumnRecord 构建列定义记录，num 由 counter 递增分配。
func BuildColumnRecord(c *counter.Counter, col *columnv1.ColumnRecord) *recordv1.Record {
	return &recordv1.Record{
		Num:       c.Inc(),
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Column{
			Column: col,
		},
	}
}

// BuildScalarRecord 构建标量数据记录。
func BuildScalarRecord(c *counter.Counter, rec *datav1.ScalarRecord) *recordv1.Record {
	return &recordv1.Record{
		Num:       c.Inc(),
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Scalar{
			Scalar: rec,
		},
	}
}

// BuildMediaRecord 构建媒体数据记录。
func BuildMediaRecord(c *counter.Counter, rec *datav1.MediaRecord) *recordv1.Record {
	return &recordv1.Record{
		Num:       c.Inc(),
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Media{
			Media: rec,
		},
	}
}

// BuildLogRecord 构建控制台日志记录：拷贝入参并分配递增的 epoch 行号，
// num 由 counter 递增分配。入参本身不被修改。
func BuildLogRecord(c *counter.Counter, epoch *counter.Counter, rec *terminalv1.LogRecord) *recordv1.Record {
	cp := proto.Clone(rec).(*terminalv1.LogRecord)
	cp.Epoch = epoch.Inc()
	return &recordv1.Record{
		Num:       c.Inc(),
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Log{
			Log: cp,
		},
	}
}

// BuildStartRecord 构建启动记录，num 固定为 -1（仅写本地 store 通道）。
func BuildStartRecord(rec *runv1.StartRecord) *recordv1.Record {
	return &recordv1.Record{
		Num:       record.NumStart,
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Start{
			Start: rec,
		},
	}
}

// BuildFinishRecord 构建结束记录，num 固定为 -2（仅写本地 store 通道）。
func BuildFinishRecord(rec *runv1.FinishRecord) *recordv1.Record {
	return &recordv1.Record{
		Num:       record.NumFinish,
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Finish{
			Finish: rec,
		},
	}
}

// BuildSaveRecord 按 rec.Type 分配编号：CUSTOM 走 counter 递增，
// 四类内部元数据使用固定槽位（Config=-3、Metadata=-4、Requirements=-5、Conda=-6，
// 可重复上报、由上传侧同号覆盖去重）。未知类型返回错误。
func BuildSaveRecord(c *counter.Counter, rec *savev1.SaveRecord) (*recordv1.Record, error) {
	var num int64
	switch rec.GetType() {
	case savev1.SaveType_SAVE_TYPE_CUSTOM:
		num = c.Inc()
	case savev1.SaveType_SAVE_TYPE_CONFIG:
		num = record.NumConfigSave
	case savev1.SaveType_SAVE_TYPE_METADATA:
		num = record.NumMetadataSave
	case savev1.SaveType_SAVE_TYPE_REQUIREMENTS:
		num = record.NumRequirementsSave
	case savev1.SaveType_SAVE_TYPE_CONDA:
		num = record.NumCondaSave
	default:
		return nil, fmt.Errorf("unsupported save type: %v", rec.GetType())
	}
	return &recordv1.Record{
		Num:       num,
		Timestamp: timestamppb.Now(),
		RecordType: &recordv1.Record_Save{
			Save: rec,
		},
	}, nil
}
