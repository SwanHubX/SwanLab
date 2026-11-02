package builder

import (
	"strings"

	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
)

// SystemKeyPrefix 是系统指标 key 的统一前缀。系统列（如硬件监控指标）
// 由此前缀标记，恢复摘要时用于区分 global_step 与 global_system_step。
const SystemKeyPrefix = "__swanlab__."

// sectionMinParts 是切分 section 所需的最少段数。
const sectionMinParts = 2

// IsSystemKey 判断指标 key 是否为系统指标 key。
func IsSystemKey(key string) bool {
	return strings.HasPrefix(key, SystemKeyPrefix)
}

// pyMod 返回非负取模结果，保证负数 sectionRule 也能正确切分。
func pyMod(a, n int) int {
	m := a % n
	if m < 0 {
		m += n
	}
	return m
}

// BuildAutoColumn 为首次出现的指标 key 自动构建列定义：
// 以 "/" 切分 key，按 sectionRule 决定归入 section 的前缀段数，
// 列类别为 CUSTOM、section 类型为 PUBLIC。
//
// 单段 key 的 section 为空；多段 key 的切分位置为
// sectionRule % (len(parts)-1) + 1。
func BuildAutoColumn(sectionRule int32, key string, typ columnv1.ColumnType) *columnv1.ColumnRecord {
	var section string
	parts := strings.Split(key, "/")
	if len(parts) >= sectionMinParts {
		cut := pyMod(int(sectionRule), len(parts)-1) + 1
		section = strings.Join(parts[:cut], "/")
	}
	return &columnv1.ColumnRecord{
		ColumnClass: columnv1.ColumnClass_COLUMN_CLASS_CUSTOM,
		ColumnKey:   key,
		ColumnType:  typ,
		SectionName: section,
		SectionType: columnv1.SectionType_SECTION_TYPE_PUBLIC,
	}
}

// BuildResumeColumn 构建 Resume 模式下从云端摘要恢复的降级列定义：
// 只恢复 key、type、class 三要素（section name 等展示信息不恢复）。
//
// media=true 时列类型用 UNSPECIFIED 占位（摘要无法得知媒体列的具体类型），
// 后续 CommitMedia 对占位列放行任意媒体类型；system 仅对标量列生效。
func BuildResumeColumn(key string, media bool, system bool) *columnv1.ColumnRecord {
	if media {
		return &columnv1.ColumnRecord{
			ColumnKey:  key,
			ColumnType: columnv1.ColumnType_COLUMN_TYPE_UNSPECIFIED,
		}
	}
	class := columnv1.ColumnClass_COLUMN_CLASS_CUSTOM
	if system {
		class = columnv1.ColumnClass_COLUMN_CLASS_SYSTEM
	}
	return &columnv1.ColumnRecord{
		ColumnKey:   key,
		ColumnType:  columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		ColumnClass: class,
	}
}

// MediumName 返回媒体类型对应的目录名，同时用作云端路径的 media/<medium>/
// 段；未知类型（含 UNSPECIFIED）返回 "unknown"。
func MediumName(typ columnv1.ColumnType) string {
	switch typ {
	case columnv1.ColumnType_COLUMN_TYPE_TEXT:
		return "text"
	case columnv1.ColumnType_COLUMN_TYPE_IMAGE:
		return "image"
	case columnv1.ColumnType_COLUMN_TYPE_AUDIO:
		return "audio"
	case columnv1.ColumnType_COLUMN_TYPE_VIDEO:
		return "video"
	case columnv1.ColumnType_COLUMN_TYPE_ECHARTS:
		return "echarts"
	case columnv1.ColumnType_COLUMN_TYPE_OBJECT3D:
		return "object3d"
	case columnv1.ColumnType_COLUMN_TYPE_MOLECULE:
		return "molecule"
	case columnv1.ColumnType_COLUMN_TYPE_HTML:
		return "html"
	default:
		return "unknown"
	}
}
