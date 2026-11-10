package manager

import (
	"context"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/budget"
	"github.com/swanhubx/swanlab/core/internal/uploader"
	operationv1 "github.com/swanhubx/swanlab/core/proto/swanlab/operation/v1"
)

func TestSessionSnapshotIncludesUploadStatus(t *testing.T) {
	up := newIdleUploader(t)
	up.Suspend()
	t.Cleanup(up.Resume)

	tracker := uploader.NewTracker()
	tracker.SetState(operationv1.CoreState_CORE_STATE_RUNNING)
	ledger := budget.New(budget.Limits{})
	if err := ledger.Acquire(context.Background(), 0, 512); err != nil {
		t.Fatalf("budget acquire: %v", err)
	}

	s := &Session{tracker: tracker, uploader: up, budget: ledger}
	stats := s.snapshot()

	if stats.GetState() != operationv1.CoreState_CORE_STATE_RUNNING {
		t.Fatalf("state = %v, want RUNNING", stats.GetState())
	}
	if !stats.GetUpload().GetPaused() {
		t.Fatal("upload.paused = false, want true")
	}
	if got := stats.GetUpload().GetBufferedBytes(); got != 512 {
		t.Fatalf("upload.buffered_bytes = %d, want 512", got)
	}
	if got := stats.GetUpload().GetPendingRecords(); got != 0 {
		t.Fatalf("upload.pending_records = %d, want 0", got)
	}
}

func TestSessionSnapshotIdleUploadStatus(t *testing.T) {
	up := newIdleUploader(t)
	tracker := uploader.NewTracker()
	s := &Session{tracker: tracker, uploader: up}
	stats := s.snapshot()

	if stats.GetUpload() == nil {
		t.Fatal("upload status should be present once the pipeline is assembled")
	}
	if stats.GetUpload().GetPaused() {
		t.Fatal("upload.paused = true without suspension")
	}
	if got := stats.GetUpload().GetBufferedBytes(); got != 0 {
		t.Fatalf("upload.buffered_bytes = %d, want 0", got)
	}
}
