package metrics

import (
	"errors"
	"io"
	"path/filepath"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/pkg/builder"
	"github.com/swanhubx/swanlab/core/internal/pkg/console"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	runv1 "github.com/swanhubx/swanlab/core/proto/swanlab/run/v1"
)

// RebuildResult 是从本地 run-*.swanlab 一次顺序扫描重建的会话状态。
type RebuildResult struct {
	Metrics          *RunMetrics
	Start            *runv1.StartRecord
	Counter          int64
	Epoch            int64
	GlobalStep       int64
	GlobalSystemStep int64
	EndSaves         []*recordv1.Record
	Records          uint64
	EndOffset        int64
	Finished         bool
}

// Rebuild 从文件头顺序扫描本地 store 文件，一次构建 RunMetrics、Counter、
// console epoch、global_step/global_system_step 与延迟上传的 SAVE_POLICY_END 记录。
//
// 遇到 ErrTruncated、ErrCorrupt 或损坏的 protobuf 载荷时在最后一条完整记录处停止，
// 告警一次并保留可信前缀；文件无法打开或文件头无效时返回错误。
func Rebuild(path string) (*RebuildResult, error) {
	r, err := store.Open(path)
	if err != nil {
		return nil, err
	}
	defer func() { _ = r.Close() }()

	res := &RebuildResult{
		Metrics:          New(),
		GlobalStep:       -1,
		GlobalSystemStep: -1,
		EndOffset:        store.HeaderLen,
	}

	base := filepath.Base(path)
	for {
		payload, scanErr := r.Scan()
		if scanErr != nil {
			switch {
			case errors.Is(scanErr, io.EOF):
				return res, nil
			case errors.Is(scanErr, store.ErrTruncated), errors.Is(scanErr, store.ErrCorrupt):
				console.Warningf("rebuild: %s stopped at offset %d after %d records: %v",
					base, res.EndOffset, res.Records, scanErr)
				return res, nil
			default:
				return nil, scanErr
			}
		}

		rec := &recordv1.Record{}
		if uerr := proto.Unmarshal(payload, rec); uerr != nil {
			console.Warningf("rebuild: %s corrupt payload at offset %d after %d records: %v",
				base, res.EndOffset, res.Records, uerr)
			return res, nil
		}

		res.Records++
		res.EndOffset = r.Position()
		applyRebuiltRecord(res, rec)
	}
}

func applyRebuiltRecord(res *RebuildResult, rec *recordv1.Record) {
	if rec.GetNum() > res.Counter {
		res.Counter = rec.GetNum()
	}

	switch record.TypeOf(rec) {
	case record.TypeStart:
		if res.Records == 1 && rec.GetNum() == record.NumStart && rec.GetStart() != nil {
			res.Start = rec.GetStart()
		}
	case record.TypeFinish:
		res.Finished = true
	case record.TypeColumn:
		applyRebuiltColumn(res.Metrics, rec.GetColumn())
	case record.TypeScalar:
		applyRebuiltScalar(res, rec)
	case record.TypeMedia:
		applyRebuiltMedia(res, rec)
	case record.TypeHistogram:
		applyRebuiltHistogram(res, rec)
	case record.TypeLog:
		if lg := rec.GetLog(); lg != nil && lg.GetEpoch() > res.Epoch {
			res.Epoch = lg.GetEpoch()
		}
	case record.TypeSave:
		if record.IsDeferredSave(rec) {
			res.EndSaves = append(res.EndSaves, rec)
		}
	}
}

func applyRebuiltColumn(m *RunMetrics, col *columnv1.ColumnRecord) {
	if col == nil || col.GetColumnKey() == "" || m.Defined(col.GetColumnKey()) {
		return
	}
	_ = m.Define(col.GetColumnKey(), col, -1, rebuildColumnKind(col.GetColumnType()))
}

func rebuildColumnKind(t columnv1.ColumnType) Kind {
	switch t {
	case columnv1.ColumnType_COLUMN_TYPE_SCALAR:
		return KindScalar
	case columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM:
		return KindHistogram
	default:
		return KindMedia
	}
}

func applyRebuiltScalar(res *RebuildResult, rec *recordv1.Record) {
	sc := rec.GetScalar()
	if sc == nil || sc.GetKey() == "" {
		return
	}
	if !res.Metrics.Defined(sc.GetKey()) {
		col := builder.BuildResumeColumn(sc.GetKey(), false, builder.IsSystemKey(sc.GetKey()))
		_ = res.Metrics.DefineScalar(sc.GetKey(), col, -1)
	}
	if accepted, err := res.Metrics.CommitScalar(sc, nil); err == nil && accepted {
		updateGlobalStep(res, sc.GetKey(), sc.GetStep(), true)
	}
}

func applyRebuiltMedia(res *RebuildResult, rec *recordv1.Record) {
	md := rec.GetMedia()
	if md == nil || md.GetKey() == "" {
		return
	}
	if !res.Metrics.Defined(md.GetKey()) {
		col := builder.BuildResumeColumn(md.GetKey(), true, false)
		_ = res.Metrics.DefineMedia(md.GetKey(), col, -1)
	}
	if accepted, err := res.Metrics.CommitMedia(md, nil); err == nil && accepted {
		updateGlobalStep(res, md.GetKey(), md.GetStep(), false)
	}
}

func applyRebuiltHistogram(res *RebuildResult, rec *recordv1.Record) {
	hg := rec.GetHistogram()
	if hg == nil || hg.GetKey() == "" {
		return
	}
	if !res.Metrics.Defined(hg.GetKey()) {
		col := builder.BuildHistogramResumeColumn(hg.GetKey())
		_ = res.Metrics.DefineHistogram(hg.GetKey(), col, -1)
	}
	if accepted, err := res.Metrics.CommitHistogram(hg, nil); err == nil && accepted {
		updateGlobalStep(res, hg.GetKey(), hg.GetStep(), false)
	}
}

func updateGlobalStep(res *RebuildResult, key string, step int64, checkSystem bool) {
	if checkSystem && builder.IsSystemKey(key) {
		if step > res.GlobalSystemStep {
			res.GlobalSystemStep = step
		}
		return
	}
	if step > res.GlobalStep {
		res.GlobalStep = step
	}
}
