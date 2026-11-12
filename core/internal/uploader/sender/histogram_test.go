package sender

import (
	"context"
	"net/http"
	"sync/atomic"
	"testing"

	"github.com/swanhubx/swanlab/core/internal/clienttest"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
)

func histogramRecord(num int64, key string, step int64, values []uint64, bins []float64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Histogram{Histogram: &datav1.HistogramRecord{
			Key:       key,
			Step:      step,
			Timestamp: testTimestamp(),
			Value:     &datav1.HistogramValue{Values: values, Bins: bins},
		}},
	}
}

func TestHistogramUploadEncodesBuckets(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewHistogram(newSenderOptions(t, cloud))

	record := histogramRecord(1, "weights", 3, []uint64{2, 5, 1}, []float64{-1, 0, 1, 2})
	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d",
			len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}

	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(requests) != 1 {
		t.Fatalf("requests = %d, want 1", len(requests))
	}
	body := decodeBody(t, requests[0].Body)
	if body["type"] != "histogram" || body["projectId"] != "proj-1" || body["experimentId"] != "exp-1" {
		t.Fatalf("envelope = %+v", body)
	}
	metrics := body["metrics"].([]any)
	if len(metrics) != 1 {
		t.Fatalf("metrics = %d, want 1", len(metrics))
	}
	metric := metrics[0].(map[string]any)
	if metric["key"] != "weights" || metric["index"] != float64(3) {
		t.Fatalf("metric = %+v", metric)
	}
	if metric["create_time"] != "2026-03-10T19:02:03.123400Z" {
		t.Fatalf("create_time = %v", metric["create_time"])
	}
	data := metric["data"].(map[string]any)
	if len(data) != 2 {
		t.Fatalf("data = %+v, want only values and bins", data)
	}
	values := data["values"].([]any)
	if len(values) != 3 || values[0] != float64(2) || values[2] != float64(1) {
		t.Fatalf("values = %v", values)
	}
	bins := data["bins"].([]any)
	if len(bins) != 4 || bins[0] != float64(-1) || bins[3] != float64(2) {
		t.Fatalf("bins = %v", bins)
	}
}

func TestHistogramUploadChunks(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewHistogram(newSenderOptions(t, cloud))

	records := make([]*recordv1.Record, 0, 130)
	for i := 0; i < 130; i++ {
		records = append(records, histogramRecord(int64(i+1), "weights", int64(i), []uint64{1}, []float64{0, 1}))
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 130 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d",
			len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(requests) != 2 {
		t.Fatalf("requests = %d, want 2 chunks", len(requests))
	}
	for i, expected := range []int{128, 2} {
		metrics := decodeBody(t, requests[i].Body)["metrics"].([]any)
		if len(metrics) != expected {
			t.Fatalf("chunk %d metrics = %d, want %d", i, len(metrics), expected)
		}
	}
}

func TestHistogramInvalidRecordsTerminate(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewHistogram(newSenderOptions(t, cloud))

	missingValue := histogramRecord(1, "weights", 0, []uint64{1}, []float64{0, 1})
	missingValue.GetHistogram().Value = nil
	missingTimestamp := histogramRecord(2, "weights", 1, []uint64{1}, []float64{0, 1})
	missingTimestamp.GetHistogram().Timestamp = nil
	badBins := histogramRecord(3, "weights", 2, []uint64{1, 2}, []float64{0, 1})
	emptyValues := histogramRecord(4, "weights", 3, nil, []float64{0})
	tooMany := make([]uint64, maxHistogramBuckets+1)
	oversized := histogramRecord(5, "weights", 4, tooMany, make([]float64, len(tooMany)+1))
	valid := histogramRecord(6, "weights", 5, []uint64{1}, []float64{0, 1})

	result := sender.Upload(context.Background(),
		[]*recordv1.Record{missingValue, missingTimestamp, badBins, emptyValues, oversized, valid})
	if len(result.Confirmed) != 6 || len(result.Remaining) != 0 {
		t.Fatalf("result = confirmed %d remaining %d", len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 2 ||
		result.Failures[0].Message != histogramBlockedMessage ||
		result.Failures[1].Message != histogramInvalidShapeMessage {
		t.Fatalf("failures = %+v", result.Failures)
	}
	if len(result.Failures[0].Records) != 2 || len(result.Failures[1].Records) != 3 {
		t.Fatalf("failure records = %d/%d, want 2/3",
			len(result.Failures[0].Records), len(result.Failures[1].Records))
	}
	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(requests) != 1 {
		t.Fatalf("requests = %d, want only the valid record", len(requests))
	}
	metrics := decodeBody(t, requests[0].Body)["metrics"].([]any)
	if len(metrics) != 1 || metrics[0].(map[string]any)["index"] != float64(5) {
		t.Fatalf("metrics = %v, want the valid record", metrics)
	}
}

func TestHistogramRetryableKeepsFailedChunk(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var calls atomic.Int32
	cloud.Handle(http.MethodPost, metricsPath, func(w http.ResponseWriter, _ *http.Request) {
		if calls.Add(1) == 1 {
			clienttest.WriteJSON(w, http.StatusOK, nil)
			return
		}
		clienttest.WriteJSON(w, http.StatusServiceUnavailable, nil)
	})
	sender := NewHistogram(newSenderOptions(t, cloud))

	records := make([]*recordv1.Record, 0, 129)
	for i := 0; i < 129; i++ {
		records = append(records, histogramRecord(int64(i+1), "weights", int64(i), []uint64{1}, []float64{0, 1}))
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 128 || len(result.Remaining) != 1 {
		t.Fatalf("result = confirmed %d remaining %d, want 128/1",
			len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 1 {
		t.Fatalf("failures = %+v", result.Failures)
	}
}

func TestHistogramDeterministicErrorContinuesNextChunk(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var calls atomic.Int32
	cloud.Handle(http.MethodPost, metricsPath, func(w http.ResponseWriter, _ *http.Request) {
		if calls.Add(1) == 1 {
			clienttest.WriteJSON(w, http.StatusUnauthorized, nil)
			return
		}
		clienttest.WriteJSON(w, http.StatusOK, nil)
	})
	sender := NewHistogram(newSenderOptions(t, cloud))

	records := make([]*recordv1.Record, 0, 129)
	for i := 0; i < 129; i++ {
		records = append(records, histogramRecord(int64(i+1), "weights", int64(i), []uint64{1}, []float64{0, 1}))
	}
	result := sender.Upload(context.Background(), records)

	// 确定性错误只终结失败的那一片，第二片必须真的发出去。
	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(requests) != 2 {
		t.Fatalf("requests = %d, want 2: the second chunk must still be attempted", len(requests))
	}
	if len(decodeBody(t, requests[1].Body)["metrics"].([]any)) != 1 {
		t.Fatalf("second chunk body does not carry the leftover record")
	}
	if len(result.Confirmed) != 129 || len(result.Remaining) != 0 {
		t.Fatalf("result = confirmed %d remaining %d, want 129/0",
			len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 1 || len(result.Failures[0].Records) != histogramChunk {
		t.Fatalf("failures = %+v, want exactly the first %d records", result.Failures, histogramChunk)
	}
}

func TestHistogramDeterministicErrorTerminates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusUnauthorized, nil)
	sender := NewHistogram(newSenderOptions(t, cloud))

	records := []*recordv1.Record{
		histogramRecord(1, "weights", 0, []uint64{1}, []float64{0, 1}),
		histogramRecord(2, "weights", 1, []uint64{1}, []float64{0, 1}),
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 2 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d",
			len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
}
