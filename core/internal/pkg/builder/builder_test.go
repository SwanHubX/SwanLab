package builder

import (
	"testing"

	"github.com/swanhubx/swanlab/core/internal/pkg/counter"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

func TestBuildColumnRecord(t *testing.T) {
	c := counter.New(0)
	col := &columnv1.ColumnRecord{ColumnKey: "loss", ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR}
	r := BuildColumnRecord(c, col)
	if r.Num != 1 {
		t.Fatalf("Num = %d, want 1", r.Num)
	}
	if r.GetColumn() != col {
		t.Fatal("Column 负载未指向入参")
	}
	if r.GetTimestamp() == nil {
		t.Fatal("Timestamp 未设置")
	}
	if got := record.TypeOf(r); got != record.TypeColumn {
		t.Fatalf("TypeOf = %q, want column", got)
	}
	if got := c.Value(); got != 1 {
		t.Fatalf("counter = %d, want 1", got)
	}
}

func TestBuildScalarAndMediaRecord(t *testing.T) {
	c := counter.New(0)
	sr := BuildScalarRecord(c, &datav1.ScalarRecord{Key: "loss", Step: 3})
	mr := BuildMediaRecord(c, &datav1.MediaRecord{Key: "img", Step: 4})
	if sr.Num != 1 || mr.Num != 2 {
		t.Fatalf("nums = %d,%d want 1,2", sr.Num, mr.Num)
	}
	if record.TypeOf(sr) != record.TypeScalar || record.TypeOf(mr) != record.TypeMedia {
		t.Fatalf("类型标签错误: %q %q", record.TypeOf(sr), record.TypeOf(mr))
	}
}

func TestBuildStartFinishRecord(t *testing.T) {
	sr := BuildStartRecord(&runv1.StartRecord{Name: "n"})
	fr := BuildFinishRecord(&runv1.FinishRecord{})
	if sr.Num != record.NumStart {
		t.Fatalf("start num = %d, want %d", sr.Num, record.NumStart)
	}
	if fr.Num != record.NumFinish {
		t.Fatalf("finish num = %d, want %d", fr.Num, record.NumFinish)
	}
	// 固定槽位不消耗 counter：Start/Finish 之后正数编号仍从 1 开始
	c := counter.New(0)
	if got := BuildColumnRecord(c, &columnv1.ColumnRecord{}).Num; got != 1 {
		t.Fatalf("start/finish 之后首个正数编号 = %d, want 1", got)
	}
}

func TestBuildLogRecord(t *testing.T) {
	c := counter.New(0)
	epoch := counter.New(0)
	in := &terminalv1.LogRecord{Line: "hello", Level: terminalv1.LogLevel_LOG_LEVEL_INFO}

	r1 := BuildLogRecord(c, epoch, in)
	r2 := BuildLogRecord(c, epoch, in)
	if r1.Num != 1 || r2.Num != 2 {
		t.Fatalf("nums = %d,%d want 1,2", r1.Num, r2.Num)
	}
	if r1.GetLog().Epoch != 1 || r2.GetLog().Epoch != 2 {
		t.Fatalf("epochs = %d,%d want 1,2", r1.GetLog().Epoch, r2.GetLog().Epoch)
	}
	// 入参未被修改（epoch 由拷贝承载）
	if in.Epoch != 0 {
		t.Fatalf("入参 LogRecord 被修改: epoch=%d", in.Epoch)
	}
	if record.TypeOf(r1) != record.TypeLog {
		t.Fatalf("TypeOf = %q, want log", record.TypeOf(r1))
	}
}

func TestBuildSaveRecord(t *testing.T) {
	c := counter.New(0)
	cases := []struct {
		typ savev1.SaveType
		num int64
	}{
		{savev1.SaveType_SAVE_TYPE_METADATA, -4},
		{savev1.SaveType_SAVE_TYPE_REQUIREMENTS, -5},
		{savev1.SaveType_SAVE_TYPE_CONDA, -6},
		{savev1.SaveType_SAVE_TYPE_CONFIG, -3},
		{savev1.SaveType_SAVE_TYPE_CUSTOM, 1}, // 首个正数编号
	}
	for i, tc := range cases {
		r, err := BuildSaveRecord(c, &savev1.SaveRecord{Name: "f", Type: tc.typ})
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		if r.Num != tc.num {
			t.Fatalf("case %d (%v): num = %d, want %d", i, tc.typ, r.Num, tc.num)
		}
	}
	// 未知类型返回错误
	if _, err := BuildSaveRecord(c, &savev1.SaveRecord{Type: savev1.SaveType(99)}); err == nil {
		t.Fatal("未知 save 类型应返回错误")
	}
}

func TestNumSequenceMirrorsPython(t *testing.T) {
	// 一次典型上报顺序：start → column → scalar → log → save → finish。
	// 断言 num 序列为 [-1, 1, 2, 3, 4, -2]（对齐 S3 冒烟的断言口径）。
	c := counter.New(0)
	epoch := counter.New(0)

	nums := make([]int64, 0, 6)
	nums = append(nums,
		BuildStartRecord(&runv1.StartRecord{}).Num,
		BuildColumnRecord(c, &columnv1.ColumnRecord{ColumnKey: "loss"}).Num,
		BuildScalarRecord(c, &datav1.ScalarRecord{Key: "loss", Step: 1}).Num,
		BuildLogRecord(c, epoch, &terminalv1.LogRecord{Line: "x"}).Num,
	)
	sr, err := BuildSaveRecord(c, &savev1.SaveRecord{Type: savev1.SaveType_SAVE_TYPE_CUSTOM})
	if err != nil {
		t.Fatal(err)
	}
	nums = append(nums, sr.Num, BuildFinishRecord(&runv1.FinishRecord{}).Num)

	want := []int64{-1, 1, 2, 3, 4, -2}
	for i := range want {
		if nums[i] != want[i] {
			t.Fatalf("num[%d] = %d, want %d (全序列 %v)", i, nums[i], want[i], nums)
		}
	}
}

func TestBuildAutoColumn(t *testing.T) {
	sc := columnv1.ColumnType_COLUMN_TYPE_SCALAR
	tests := []struct {
		name       string
		rule       int32
		key        string
		wantSecond string
	}{
		// cut = rule % (len(parts)-1) + 1
		{"单段无 section", 0, "loss", ""},
		{"两段 rule=0", 0, "train/loss", "train"},
		{"两段 rule=1", 1, "train/loss", "train"},
		{"三段 rule=0", 0, "a/b/c", "a"},
		{"三段 rule=1", 1, "a/b/c", "a/b"},
		{"三段 rule=2", 2, "a/b/c", "a"},
		{"三段 rule=3 回绕", 3, "a/b/c", "a/b"},
		{"四段 rule=2", 2, "x/y/z/w", "x/y/z"},
		{"尾斜杠", 1, "train/", "train"},
		{"负数 rule 取模", -1, "a/b/c", "a/b"},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			col := BuildAutoColumn(tc.rule, tc.key, sc)
			if col.SectionName != tc.wantSecond {
				t.Fatalf("SectionName = %q, want %q", col.SectionName, tc.wantSecond)
			}
			if col.ColumnKey != tc.key {
				t.Fatalf("ColumnKey = %q, want %q", col.ColumnKey, tc.key)
			}
			if col.ColumnType != sc {
				t.Fatalf("ColumnType = %v", col.ColumnType)
			}
			if col.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_CUSTOM {
				t.Fatalf("ColumnClass = %v", col.ColumnClass)
			}
			if col.SectionType != columnv1.SectionType_SECTION_TYPE_PUBLIC {
				t.Fatalf("SectionType = %v", col.SectionType)
			}
		})
	}
}

func TestBuildResumeColumn(t *testing.T) {
	media := BuildResumeColumn("pic", true, false)
	if media.ColumnType != columnv1.ColumnType_COLUMN_TYPE_UNSPECIFIED {
		t.Fatalf("media resume 列类型 = %v, want UNSPECIFIED", media.ColumnType)
	}
	if media.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_UNSPECIFIED {
		t.Fatalf("media resume 列类别 = %v", media.ColumnClass)
	}

	custom := BuildResumeColumn("loss", false, false)
	if custom.ColumnType != columnv1.ColumnType_COLUMN_TYPE_SCALAR {
		t.Fatalf("scalar resume 列类型 = %v", custom.ColumnType)
	}
	if custom.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_CUSTOM {
		t.Fatalf("custom resume 列类别 = %v", custom.ColumnClass)
	}

	sys := BuildResumeColumn("__swanlab__.gpu", false, true)
	if sys.ColumnClass != columnv1.ColumnClass_COLUMN_CLASS_SYSTEM {
		t.Fatalf("system resume 列类别 = %v", sys.ColumnClass)
	}
	if !IsSystemKey(sys.ColumnKey) {
		t.Fatal("IsSystemKey 未识别系统前缀")
	}
}

func TestMediumName(t *testing.T) {
	tests := []struct {
		typ  columnv1.ColumnType
		want string
	}{
		{columnv1.ColumnType_COLUMN_TYPE_TEXT, "text"},
		{columnv1.ColumnType_COLUMN_TYPE_IMAGE, "image"},
		{columnv1.ColumnType_COLUMN_TYPE_AUDIO, "audio"},
		{columnv1.ColumnType_COLUMN_TYPE_VIDEO, "video"},
		{columnv1.ColumnType_COLUMN_TYPE_ECHARTS, "echarts"},
		{columnv1.ColumnType_COLUMN_TYPE_OBJECT3D, "object3d"},
		{columnv1.ColumnType_COLUMN_TYPE_MOLECULE, "molecule"},
		{columnv1.ColumnType_COLUMN_TYPE_HTML, "html"},
		{columnv1.ColumnType_COLUMN_TYPE_UNSPECIFIED, "unknown"},
		{columnv1.ColumnType_COLUMN_TYPE_SCALAR, "unknown"},
	}
	for _, tc := range tests {
		if got := MediumName(tc.typ); got != tc.want {
			t.Fatalf("MediumName(%v) = %q, want %q", tc.typ, got, tc.want)
		}
	}
}

func TestGenerateNameAndColor(t *testing.T) {
	names := map[int]string{
		0:   "quick-swan-0",
		5:   "old-dragon-5",
		19:  "fresh-zebra-19",
		20:  "quick-swan-20",
		128: "happy-goat-128",
		41:  "slow-rat-41",
	}
	for slug, want := range names {
		if got := GenerateName(slug); got != want {
			t.Fatalf("GenerateName(%d) = %q, want %q", slug, got, want)
		}
	}
	colors := map[int]string{
		0:  "#528d59",
		15: "#989fa3",
		16: "#528d59", // 回绕
		21: "#dfb142",
	}
	for slug, want := range colors {
		if got := GenerateColor(slug); got != want {
			t.Fatalf("GenerateColor(%d) = %q, want %q", slug, got, want)
		}
	}
	// 负数 slug 取模（-1 % 20 = 19）
	if got := GenerateName(-1); got != "fresh-zebra--1" {
		t.Fatalf("GenerateName(-1) = %q", got)
	}
}

func TestRecordsShareCounterMonotonic(t *testing.T) {
	// 同一 counter 下混用各类构建器，正数编号必须连续无洞
	c := counter.New(0)
	var nums []int64
	for i := 0; i < 100; i++ {
		switch i % 4 {
		case 0:
			nums = append(nums, BuildColumnRecord(c, &columnv1.ColumnRecord{}).Num)
		case 1:
			nums = append(nums, BuildScalarRecord(c, &datav1.ScalarRecord{}).Num)
		case 2:
			nums = append(nums, BuildMediaRecord(c, &datav1.MediaRecord{}).Num)
		case 3:
			r, err := BuildSaveRecord(c, &savev1.SaveRecord{})
			if err != nil {
				t.Fatal(err)
			}
			nums = append(nums, r.Num)
		}
	}
	for i, n := range nums {
		if n != int64(i+1) {
			t.Fatalf("nums[%d] = %d, 编号不连续: %v", i, n, nums)
		}
	}
}
