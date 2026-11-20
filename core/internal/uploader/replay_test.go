package uploader

import (
	"errors"
	"io"
	"path/filepath"
	"testing"
	"time"

	"google.golang.org/protobuf/proto"

	"github.com/swanhubx/swanlab/core/internal/pkg/record"
	"github.com/swanhubx/swanlab/core/internal/store"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

// writeStore 把记录按 run 期口径序列化落盘。
func writeStore(t *testing.T, path string, records []*recordv1.Record) {
	t.Helper()
	writer, err := store.Create(path)
	if err != nil {
		t.Fatalf("create store: %v", err)
	}
	for _, rec := range records {
		payload, err := proto.Marshal(rec)
		if err != nil {
			t.Fatalf("marshal record: %v", err)
		}
		if _, err := writer.WriteRecord(payload); err != nil {
			t.Fatalf("write record: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("close store: %v", err)
	}
}

// readStore 顺序读回整份 store，与 T4 回填的读取口径一致。
func readStore(t *testing.T, path string) []*recordv1.Record {
	t.Helper()
	reader, err := store.Open(path)
	if err != nil {
		t.Fatalf("open store: %v", err)
	}
	defer func() { _ = reader.Close() }()

	var out []*recordv1.Record
	for {
		payload, err := reader.Scan()
		if errors.Is(err, io.EOF) {
			return out
		}
		if err != nil {
			t.Fatalf("scan store: %v", err)
		}
		rec := &recordv1.Record{}
		if err := proto.Unmarshal(payload, rec); err != nil {
			t.Fatalf("unmarshal record: %v", err)
		}
		out = append(out, rec)
	}
}

// TestStoreReplayRoutesHistogram 固定 T4 回填依赖的路由契约：从本地 store 读回的
// 记录按 TypeOf + dispatchOrder 分发，histogram 必须进 histogram sender，
// 不能被当作未知类型终结。
func TestStoreReplayRoutesHistogram(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run-replay.swanlab")
	written := []*recordv1.Record{
		newColumn(1, "loss"),
		newScalar(2, "loss"),
		newHistogram(3, "weights"),
		newHistogram(4, "weights"),
		newScalar(5, "acc"),
	}
	writeStore(t, path, written)

	replayed := readStore(t, path)
	if len(replayed) != len(written) {
		t.Fatalf("replayed = %d records, want %d", len(replayed), len(written))
	}
	for i, rec := range replayed {
		if record.TypeOf(rec) != record.TypeOf(written[i]) || rec.GetNum() != written[i].GetNum() {
			t.Fatalf("record %d = num %d type %s, want num %d type %s",
				i, rec.GetNum(), record.TypeOf(rec), written[i].GetNum(), record.TypeOf(written[i]))
		}
	}

	columnCalls, columnSender := recordingSender(4)
	scalarCalls, scalarSender := recordingSender(4)
	histCalls, histSender := recordingSender(4)
	f := newUplFixture(t, map[record.Type]Sender{
		record.TypeColumn:    columnSender,
		record.TypeScalar:    scalarSender,
		record.TypeHistogram: histSender,
	}, 8, time.Millisecond)
	f.feed(t, replayed...)
	close(f.in)
	waitDone(t, f.u)

	if batch := recvBatch(t, columnCalls); len(batch) != 1 || record.TypeOf(batch[0]) != record.TypeColumn {
		t.Fatalf("column batch = %v", batch)
	}
	if batch := recvBatch(t, histCalls); len(batch) != 2 ||
		batch[0].GetNum() != 3 || batch[1].GetNum() != 4 {
		t.Fatalf("histogram batch = %v, want the two replayed histogram records in batch order", batch)
	}
	if batch := recvBatch(t, scalarCalls); len(batch) != 2 {
		t.Fatalf("scalar batch = %v", batch)
	}
	if got := f.tracker.Snapshot().GetUploadedRecords(); got != int64(len(written)) {
		t.Fatalf("uploaded = %d, want %d: replayed records must not be dropped", got, len(written))
	}
}

// TestStoreReplayWithoutHistogramSenderTerminates 保证未注册 histogram sender 时记录被
// 显式终结并计入进度，而不是留在 pending 里永远不上传。
func TestStoreReplayWithoutHistogramSenderTerminates(t *testing.T) {
	path := filepath.Join(t.TempDir(), "run-no-sender.swanlab")
	writeStore(t, path, []*recordv1.Record{newHistogram(1, "weights")})

	f := newUplFixture(t, map[record.Type]Sender{
		record.TypeScalar: senderFunc(successUpload),
	}, 4, time.Millisecond)
	f.feed(t, readStore(t, path)...)
	close(f.in)
	waitDone(t, f.u)

	if got := f.tracker.Snapshot().GetUploadedRecords(); got != 1 {
		t.Fatalf("uploaded = %d, want 1: unsupported types must still be terminated", got)
	}
}
