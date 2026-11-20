package buffer

import (
	"context"
	"errors"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/budget"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// pausedIntake 构造暂停于 seq 3 的会话：前两条已上传并确认，第三条只落本地。
func pausedIntake(t *testing.T) (*Intake, *budget.Budget) {
	t.Helper()
	sink := newStoreSink(16)
	t.Cleanup(sink.stop)

	b := budget.New(budget.Limits{Records: 2})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 8})
	mustCommit(t, in, scalarRecord(1))
	mustCommit(t, in, scalarRecord(2))
	mustCommit(t, in, scalarRecord(3))
	drainUpload(in)

	records, bytes := b.Usage()
	b.Release(records, bytes)
	return in, b
}

func TestHandoverSendsTailAndResumes(t *testing.T) {
	in, b := pausedIntake(t)

	var synced bool
	err := in.Handover(
		context.Background(),
		func(context.Context) error {
			synced = true
			return nil
		},
		func() ([]*recordv1.Record, uint64, error) {
			if !synced {
				t.Error("fill must run after sync")
			}
			return []*recordv1.Record{scalarRecord(3)}, 3, nil
		},
	)
	if err != nil {
		t.Fatalf("Handover: %v", err)
	}
	if paused, seq := in.PauseState(); paused || seq != 0 {
		t.Fatalf("PauseState = (%v, %d), want resumed", paused, seq)
	}
	if got := recordNums(drainUpload(in)); len(got) != 1 || got[0] != 3 {
		t.Fatalf("upload = %v, want record 3", got)
	}

	// 恢复后提交的记录直接上传，不需要再次补投。
	mustCommit(t, in, scalarRecord(4))
	if got := recordNums(drainUpload(in)); len(got) != 1 || got[0] != 4 {
		t.Fatalf("upload after resume = %v, want record 4", got)
	}
	if records, _ := b.Usage(); records != 2 {
		t.Fatalf("budget usage = %d, want 2", records)
	}
}

func TestHandoverStaleTailKeepsPause(t *testing.T) {
	in, _ := pausedIntake(t)

	err := in.Handover(context.Background(), nil, func() ([]*recordv1.Record, uint64, error) {
		return []*recordv1.Record{scalarRecord(3)}, 2, nil
	})
	if !errors.Is(err, ErrBackfillStale) {
		t.Fatalf("Handover err = %v, want ErrBackfillStale", err)
	}
	if paused, seq := in.PauseState(); !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload = %d records, want 0 on stale handover", got)
	}
}

func TestHandoverSyncFailureKeepsPause(t *testing.T) {
	in, _ := pausedIntake(t)
	flushErr := errors.New("flush failed")
	filled := false

	err := in.Handover(
		context.Background(),
		func(context.Context) error { return flushErr },
		func() ([]*recordv1.Record, uint64, error) {
			filled = true
			return nil, 3, nil
		},
	)
	if !errors.Is(err, flushErr) {
		t.Fatalf("Handover err = %v, want flush error", err)
	}
	if filled {
		t.Fatal("fill must not run when sync fails")
	}
	if paused, seq := in.PauseState(); !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}
}

func TestHandoverBlockedWithoutUploadCapacity(t *testing.T) {
	sink := newStoreSink(16)
	defer sink.stop()

	b := budget.New(budget.Limits{Records: 1})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 2})
	mustCommit(t, in, scalarRecord(1))
	mustCommit(t, in, scalarRecord(2)) // 暂停于 seq 2
	mustCommit(t, in, scalarRecord(3))
	mustCommit(t, in, scalarRecord(4))
	drainUpload(in)
	records, bytes := b.Usage()
	b.Release(records, bytes)

	tail := []*recordv1.Record{scalarRecord(2), scalarRecord(3), scalarRecord(4)}
	err := in.Handover(context.Background(), nil, func() ([]*recordv1.Record, uint64, error) {
		return tail, 4, nil
	})
	if !errors.Is(err, ErrBackfillBlocked) {
		t.Fatalf("Handover err = %v, want ErrBackfillBlocked", err)
	}
	if paused, seq := in.PauseState(); !paused || seq != 2 {
		t.Fatalf("PauseState = (%v, %d), want (true, 2)", paused, seq)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload = %d records, want 0 when capacity is short", got)
	}
	if got, _ := b.Usage(); got != 0 {
		t.Fatalf("budget usage = %d, want 0 when handover is blocked", got)
	}
}

func TestHandoverAfterClose(t *testing.T) {
	in, _ := pausedIntake(t)
	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}
	err := in.Handover(context.Background(), nil, func() ([]*recordv1.Record, uint64, error) {
		return nil, 3, nil
	})
	if !errors.Is(err, ErrClosed) {
		t.Fatalf("Handover err = %v, want ErrClosed", err)
	}
}

// 交接持锁期间提交的记录必须排在尾部之后，不能被重复投递或插队。
func TestHandoverSerializesConcurrentCommit(t *testing.T) {
	in, b := pausedIntake(t)

	entered := make(chan struct{})
	release := make(chan struct{})
	handoverDone := make(chan error, 1)
	go func() {
		handoverDone <- in.Handover(
			context.Background(),
			func(context.Context) error {
				close(entered)
				<-release
				return nil
			},
			func() ([]*recordv1.Record, uint64, error) {
				return []*recordv1.Record{scalarRecord(3)}, 3, nil
			},
		)
	}()
	<-entered

	commitDone := make(chan error, 1)
	go func() {
		rec := scalarRecord(4)
		tk, err := in.Reserve(context.Background(), demandFor([]*recordv1.Record{rec}))
		if err != nil {
			commitDone <- err
			return
		}
		commitDone <- in.Commit(tk, []*recordv1.Record{rec})
	}()

	close(release)
	if err := <-handoverDone; err != nil {
		t.Fatalf("Handover: %v", err)
	}
	if err := <-commitDone; err != nil {
		t.Fatalf("Commit: %v", err)
	}
	if got := recordNums(drainUpload(in)); len(got) != 2 || got[0] != 3 || got[1] != 4 {
		t.Fatalf("upload order = %v, want 3,4", got)
	}
	if records, _ := b.Usage(); records != 2 {
		t.Fatalf("budget usage = %d, want 2", records)
	}
}
