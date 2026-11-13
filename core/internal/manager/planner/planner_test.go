package planner

import (
	"errors"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/config"
	"github.com/swanhubx/swanlab/core/internal/metrics"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	savev1 "github.com/swanhubx/swanlab/core/proto/swanlab/save/v1"
)

func TestDemandFor(t *testing.T) {
	scalarRec := &recordv1.Record{
		Num: 1,
		RecordType: &recordv1.Record_Scalar{
			Scalar: &datav1.ScalarRecord{Key: "loss", Step: 0},
		},
	}
	deferredRec := &recordv1.Record{
		Num: 2,
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{
				Type:   savev1.SaveType_SAVE_TYPE_CUSTOM,
				Policy: savev1.SavePolicy_SAVE_POLICY_END,
			},
		},
	}
	liveSaveRec := &recordv1.Record{
		Num: 3,
		RecordType: &recordv1.Record_Save{
			Save: &savev1.SaveRecord{
				Type:   savev1.SaveType_SAVE_TYPE_CUSTOM,
				Policy: savev1.SavePolicy_SAVE_POLICY_LIVE,
			},
		},
	}
	records := []*recordv1.Record{scalarRec, deferredRec, liveSaveRec}

	online := DemandFor(records, config.Settings{Mode: config.ModeOnline})
	wantBytes := record.Size(scalarRec) + record.Size(liveSaveRec)
	if online.Store != 3 || online.Upload != 2 || online.Deferred != 1 || online.Bytes != wantBytes {
		t.Fatalf("online DemandFor = %+v, want Store=3 Upload=2 Deferred=1 Bytes=%d", online, wantBytes)
	}

	skipStore := DemandFor(records, config.Settings{Mode: config.ModeOnline, SkipStore: true})
	if skipStore.Store != 0 || skipStore.Upload != 2 || skipStore.Deferred != 1 || skipStore.Bytes != wantBytes {
		t.Fatalf("skipStore DemandFor = %+v, want Store=0 Upload=2 Deferred=1 Bytes=%d", skipStore, wantBytes)
	}

	offline := DemandFor(records, config.Settings{Mode: config.ModeOffline})
	if offline.Store != 3 || offline.Upload != 0 || offline.Deferred != 0 || offline.Bytes != 0 {
		t.Fatalf("offline DemandFor = %+v, want Store=3 Upload=0 Deferred=0 Bytes=0", offline)
	}
}

func TestMetricPlannerDeduplicatesBatchAndChecksTypes(t *testing.T) {
	m := metrics.New()
	p := NewMetricPlanner(m, 0)

	first := p.PrecheckScalar(&datav1.ScalarRecord{
		Key:  "loss",
		Type: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		Step: 1,
	})
	if !first.NeedsColumn || first.Decision != metrics.DecisionAccept || first.Err != nil {
		t.Fatalf("first PrecheckScalar = %+v, want NeedsColumn=true DecisionAccept", first)
	}
	if p.Column("loss") == nil {
		t.Fatal("Column(loss) = nil, want auto-generated column")
	}

	dup := p.PrecheckScalar(&datav1.ScalarRecord{
		Key:  "loss",
		Type: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		Step: 1,
	})
	if dup.NeedsColumn || dup.Decision != metrics.DecisionSkip || dup.Err != nil {
		t.Fatalf("duplicate step PrecheckScalar = %+v, want NeedsColumn=false DecisionSkip", dup)
	}

	mismatch := p.PrecheckMedia(&datav1.MediaRecord{
		Key:  "loss",
		Type: columnv1.ColumnType_COLUMN_TYPE_IMAGE,
		Step: 2,
	})
	if !errors.Is(mismatch.Err, metrics.ErrTypeMismatch) {
		t.Fatalf("mismatch Err = %v, want metrics.ErrTypeMismatch", mismatch.Err)
	}

	// Reserve 之前 metrics 保持只读，DefineAll 之后才落入状态机。
	if m.Defined("loss") {
		t.Fatal("metrics.Defined(loss) = true before DefineAll")
	}
	if err := p.DefineAll(); err != nil {
		t.Fatalf("DefineAll: %v", err)
	}
	if !m.Defined("loss") {
		t.Fatal("metrics.Defined(loss) = false after DefineAll")
	}
}

func TestMetricPlannerPendingStepBounds(t *testing.T) {
	p := NewMetricPlanner(metrics.New(), 0)
	first := p.PrecheckHistogram(&datav1.HistogramRecord{
		Key:  "hist",
		Type: columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM,
		Step: 0,
	})
	if !first.NeedsColumn || first.Decision != metrics.DecisionAccept || first.Err != nil {
		t.Fatalf("first PrecheckHistogram = %+v", first)
	}

	negOne := p.PrecheckHistogram(&datav1.HistogramRecord{
		Key:  "hist",
		Type: columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM,
		Step: -1,
	})
	if negOne.Decision != metrics.DecisionSkip || negOne.Err != nil {
		t.Fatalf("step -1 = %+v, want DecisionSkip without error", negOne)
	}

	negTwo := p.PrecheckHistogram(&datav1.HistogramRecord{
		Key:  "hist",
		Type: columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM,
		Step: -2,
	})
	if negTwo.Decision != metrics.DecisionSkip || negTwo.Err != nil {
		t.Fatalf("step -2 = %+v, want DecisionSkip without error", negTwo)
	}

	overLimit := p.PrecheckHistogram(&datav1.HistogramRecord{
		Key:  "hist",
		Type: columnv1.ColumnType_COLUMN_TYPE_HISTOGRAM,
		Step: metrics.MaxStep + 1,
	})
	if !errors.Is(overLimit.Err, metrics.ErrStepLimit) {
		t.Fatalf("overLimit Err = %v, want ErrStepLimit", overLimit.Err)
	}
}
