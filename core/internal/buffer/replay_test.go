package buffer

import (
	"context"
	"errors"
	"sync/atomic"
	"testing"
	"time"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

func TestReplaySendsPausedRecordsInOrder(t *testing.T) {
	sink := newStoreSink(16)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 2})
	var admitted atomic.Int64
	in := newIntake(t, Options{
		Store:       sink.ch,
		Budget:      b,
		UploadQueue: 8,
		OnAdmit:     func(n int64) { admitted.Add(n) },
	})

	mustCommit(t, in, scalarRecord(1))
	mustCommit(t, in, scalarRecord(2))
	mustCommit(t, in, scalarRecord(3)) // 预算满，seq 3 起只落本地
	if paused, seq := in.PauseState(); !paused || seq != 3 {
		t.Fatalf("PauseState = (%v, %d), want (true, 3)", paused, seq)
	}
	if got := len(drainUpload(in)); got != 2 {
		t.Fatalf("upload before replay = %d, want 2", got)
	}
	admittedBefore := admitted.Load()

	// 模拟前两条记录上传确认，归还配额。
	records, bytes := b.Usage()
	b.Release(records, bytes)

	replayed := []*recordv1.Record{scalarRecord(3), scalarRecord(4)}
	sent, err := in.Replay(context.Background(), replayed)
	if err != nil {
		t.Fatalf("Replay: %v", err)
	}
	if sent != len(replayed) {
		t.Fatalf("replayed %d records, want %d", sent, len(replayed))
	}
	if got := admitted.Load(); got != admittedBefore {
		t.Fatalf("admitted = %d, want %d: replay must not count records twice", got, admittedBefore)
	}
	if used, _ := b.Usage(); used != 2 {
		t.Fatalf("budget usage = %d, want 2", used)
	}

	// 补投记录先于恢复后的实时记录进入上传通道。
	records, bytes = b.Usage()
	b.Release(records, bytes)
	in.ResumeUpload()
	mustCommit(t, in, scalarRecord(5))

	got := drainUpload(in)
	if len(got) != 3 || got[0].GetNum() != 3 || got[1].GetNum() != 4 || got[2].GetNum() != 5 {
		t.Fatalf("upload order = %v, want 3,4,5", recordNums(got))
	}
}

func TestReplayRequiresPause(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	mustCommit(t, in, scalarRecord(1))
	drainUpload(in)
	records, _ := b.Usage()

	if _, err := in.Replay(context.Background(), []*recordv1.Record{scalarRecord(2)}); !errors.Is(err, ErrNotPaused) {
		t.Fatalf("Replay without pause err = %v, want ErrNotPaused", err)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload = %d records, want 0", got)
	}
	if got, _ := b.Usage(); got != records {
		t.Fatalf("budget usage = %d, want %d", got, records)
	}

	if err := in.Close(context.Background()); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if _, err := in.Replay(context.Background(), []*recordv1.Record{scalarRecord(3)}); !errors.Is(err, ErrClosed) {
		t.Fatalf("Replay after close err = %v, want ErrClosed", err)
	}
}

func TestReplayWaitsForBudget(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 1})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	mustCommit(t, in, scalarRecord(1)) // 已上传但尚未确认，占满预算
	mustCommit(t, in, scalarRecord(2)) // 暂停于 seq 2
	drainUpload(in)

	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
	defer cancel()
	sent, err := in.Replay(ctx, []*recordv1.Record{scalarRecord(3)})
	if !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("Replay err = %v, want DeadlineExceeded", err)
	}
	if sent != 0 {
		t.Fatalf("replayed %d records while budget was full, want 0", sent)
	}
	if got := len(drainUpload(in)); got != 0 {
		t.Fatalf("upload = %d records, want 0 while budget is full", got)
	}
	if got, _ := b.Usage(); got != 1 {
		t.Fatalf("budget usage = %d, want 1", got)
	}

	records, bytes := b.Usage()
	b.Release(records, bytes)
	if _, err := in.Replay(context.Background(), []*recordv1.Record{scalarRecord(3)}); err != nil {
		t.Fatalf("Replay after release: %v", err)
	}
	if got := drainUpload(in); len(got) != 1 || got[0].GetNum() != 3 {
		t.Fatalf("upload = %v, want record 3", recordNums(got))
	}
}

// 单条记录的字节数超过预算上限时不能无限等待：按回收补投的口径超额投递。
func TestReplayOversizedRecordExceedsByteBudget(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Bytes: 1})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 4})

	mustCommit(t, in, scalarRecord(1))
	if paused, seq := in.PauseState(); !paused || seq != 1 {
		t.Fatalf("PauseState = (%v, %d), want (true, 1)", paused, seq)
	}
	if _, err := in.Replay(context.Background(), []*recordv1.Record{scalarRecord(2)}); err != nil {
		t.Fatalf("Replay: %v", err)
	}
	got := drainUpload(in)
	if len(got) != 1 || got[0].GetNum() != 2 {
		t.Fatalf("upload = %v, want record 2", recordNums(got))
	}
	if _, bytes := b.Usage(); bytes <= 1 {
		t.Fatalf("budget bytes = %d, want the oversized record charged", bytes)
	}
}

func TestOnPauseNotifiesEachInterval(t *testing.T) {
	sink := newStoreSink(16)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 1})
	pauses := make(chan uint64, 4)
	in := newIntake(t, Options{
		Store:       sink.ch,
		Budget:      b,
		UploadQueue: 8,
		OnPause:     func(seq uint64) { pauses <- seq },
	})

	mustCommit(t, in, scalarRecord(1)) // 上传
	mustCommit(t, in, scalarRecord(2)) // 暂停于 seq 2
	mustCommit(t, in, scalarRecord(3)) // 已在暂停，不重复通知
	drainUpload(in)

	records, bytes := b.Usage()
	b.Release(records, bytes)
	in.ResumeUpload()
	mustCommit(t, in, scalarRecord(4)) // 上传
	mustCommit(t, in, scalarRecord(5)) // 新暂停区间，seq 5
	drainUpload(in)

	if got := <-pauses; got != 2 {
		t.Fatalf("first pause = %d, want 2", got)
	}
	if got := <-pauses; got != 5 {
		t.Fatalf("second pause = %d, want 5", got)
	}
	select {
	case got := <-pauses:
		t.Fatalf("unexpected extra pause notification = %d", got)
	default:
	}
}

// 预算只够一部分记录时，先投递可用部分并返回条数，剩余记录由调用方重试。
func TestReplayStopsAtBudgetBoundary(t *testing.T) {
	sink := newStoreSink(8)
	defer sink.stop()
	b := budget.New(budget.Limits{Records: 2})
	in := newIntake(t, Options{Store: sink.ch, Budget: b, UploadQueue: 8})

	mustCommit(t, in, scalarRecord(1))
	mustCommit(t, in, scalarRecord(2))
	mustCommit(t, in, scalarRecord(3)) // 暂停于 seq 3
	drainUpload(in)
	// 模拟第一条记录上传确认，留出一条记录的预算。
	b.Release(1, record.Size(scalarRecord(1)))

	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Millisecond)
	defer cancel()
	sent, err := in.Replay(ctx, []*recordv1.Record{scalarRecord(4), scalarRecord(5), scalarRecord(6)})
	if sent != 1 {
		t.Fatalf("replayed %d records, want 1 before the budget ran out", sent)
	}
	if !errors.Is(err, context.DeadlineExceeded) {
		t.Fatalf("Replay err = %v, want DeadlineExceeded", err)
	}
	got := drainUpload(in)
	if len(got) != 1 || got[0].GetNum() != 4 {
		t.Fatalf("upload = %v, want record 4 only", recordNums(got))
	}
}
