package sender

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"math"
	"net/http"
	"strings"
	"sync/atomic"
	"testing"

	"google.golang.org/protobuf/types/known/timestamppb"

	"github.com/swanhubx/swanlab/core/internal/api"
	"github.com/swanhubx/swanlab/core/internal/client"
	"github.com/swanhubx/swanlab/core/internal/clienttest"
	columnv1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/column/v1"
	datav1 "github.com/swanhubx/swanlab/core/proto/swanlab/metric/data/v1"
	recordv1 "github.com/swanhubx/swanlab/core/proto/swanlab/record/v1"
	terminalv1 "github.com/swanhubx/swanlab/core/proto/swanlab/terminal/v1"
)

const (
	verifyPath  = "/api/auth/verify"
	seriesPath  = "/api/projects/alice/demo/series"
	metricsPath = "/api/house/metrics"
)

// newTestAPI 启动假后端并构造带鉴权的 api 客户端。
func newTestAPI(t *testing.T, cloud *clienttest.Cloud) *api.Client {
	t.Helper()
	cloud.ServeJSON(http.MethodGet, verifyPath, http.StatusOK, map[string]any{
		"uid":       42,
		"username":  "alice",
		"createdAt": "2026-01-01T00:00:00.000Z",
	})
	transport, err := client.New(context.Background(), cloud.URL(), "smoke-key", client.Options{})
	if err != nil {
		t.Fatalf("create client: %v", err)
	}
	return api.New(transport)
}

func newSenderOptions(t *testing.T, cloud *clienttest.Cloud) Options {
	return Options{
		API:          newTestAPI(t, cloud),
		Username:     "alice",
		Project:      "demo",
		ProjectID:    "proj-1",
		ExperimentID: "exp-1",
	}
}

func testTimestamp() *timestamppb.Timestamp { return timestamppb.New(clienttest.FixedNow) }

func columnRecord(num int64, key string) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Column{Column: &columnv1.ColumnRecord{
			ColumnKey:  key,
			ColumnType: columnv1.ColumnType_COLUMN_TYPE_SCALAR,
		}},
	}
}

func scalarRecord(num int64, key string, step int64, value float64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Scalar{Scalar: &datav1.ScalarRecord{
			Key:       key,
			Step:      step,
			Timestamp: testTimestamp(),
			Value:     &datav1.ScalarValue{Number: value},
		}},
	}
}

func logRecord(num int64, line string, level terminalv1.LogLevel, epoch int64) *recordv1.Record {
	return &recordv1.Record{
		Num: num,
		RecordType: &recordv1.Record_Log{Log: &terminalv1.LogRecord{
			Line:      line,
			Level:     level,
			Epoch:     epoch,
			Timestamp: testTimestamp(),
		}},
	}
}

func decodeBody(t *testing.T, raw []byte) map[string]any {
	t.Helper()
	var body map[string]any
	if err := json.Unmarshal(raw, &body); err != nil {
		t.Fatalf("decode body %s: %v", raw, err)
	}
	return body
}

func TestColumnUploadChunks(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, seriesPath, http.StatusOK, nil)
	sender := NewColumn(newSenderOptions(t, cloud))

	records := make([]*recordv1.Record, 0, 3001)
	for i := 0; i < 3001; i++ {
		records = append(records, columnRecord(int64(i+1), "k"))
	}
	result := sender.Upload(context.Background(), records)

	if len(result.Confirmed) != 3001 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
	requests := cloud.RequestsFor(http.MethodPost, seriesPath)
	if len(requests) != 2 {
		t.Fatalf("requests = %d, want 2 chunks", len(requests))
	}
	for i, expected := range []int{3000, 1} {
		body := decodeBody(t, requests[i].Body)
		series, _ := body["series"].([]any)
		if len(series) != expected {
			t.Fatalf("chunk[%d] series = %d, want %d", i, len(series), expected)
		}
	}
}

func TestColumnEncodeFailureSkipsHTTP(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, seriesPath, http.StatusOK, nil)
	sender := NewColumn(newSenderOptions(t, cloud))

	valid := columnRecord(1, "ok")
	invalid := columnRecord(2, "bad")
	invalid.GetColumn().ColumnType = columnv1.ColumnType_COLUMN_TYPE_UNSPECIFIED

	result := sender.Upload(context.Background(), []*recordv1.Record{valid, invalid})
	if len(result.Confirmed) != 2 || len(result.Remaining) != 0 {
		t.Fatalf("result = confirmed %d remaining %d, want all confirmed", len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 1 || !strings.HasPrefix(result.Failures[0].Message, "failed to encode column") {
		t.Fatalf("failures = %+v", result.Failures)
	}
	if got := cloud.CountFor(http.MethodPost, seriesPath); got != 0 {
		t.Fatalf("series requests = %d, want 0", got)
	}
}

func TestColumnRetryableStopsAtFailedChunk(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var calls atomic.Int32
	cloud.Handle(http.MethodPost, seriesPath, func(w http.ResponseWriter, _ *http.Request) {
		if calls.Add(1) == 1 {
			clienttest.WriteJSON(w, http.StatusOK, nil)
			return
		}
		clienttest.WriteJSON(w, http.StatusServiceUnavailable, nil)
	})
	sender := NewColumn(newSenderOptions(t, cloud))

	records := make([]*recordv1.Record, 0, 3001)
	for i := 0; i < 3001; i++ {
		records = append(records, columnRecord(int64(i+1), "k"))
	}
	result := sender.Upload(context.Background(), records)

	if len(result.Confirmed) != 3000 || len(result.Remaining) != 1 {
		t.Fatalf("result = confirmed %d remaining %d, want 3000/1", len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 1 || result.Failures[0].Records[0] != records[3000] {
		t.Fatalf("failures = %+v, want the second chunk", result.Failures)
	}
}

func TestColumnNonRetryableTerminatesRest(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	var calls atomic.Int32
	cloud.Handle(http.MethodPost, seriesPath, func(w http.ResponseWriter, _ *http.Request) {
		if calls.Add(1) == 1 {
			clienttest.WriteJSON(w, http.StatusOK, nil)
			return
		}
		clienttest.WriteJSON(w, http.StatusBadRequest, nil)
	})
	sender := NewColumn(newSenderOptions(t, cloud))

	records := make([]*recordv1.Record, 0, 3001)
	for i := 0; i < 3001; i++ {
		records = append(records, columnRecord(int64(i+1), "k"))
	}
	result := sender.Upload(context.Background(), records)

	if len(result.Confirmed) != 3001 || len(result.Remaining) != 0 {
		t.Fatalf("result = confirmed %d remaining %d, want all confirmed", len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 1 {
		t.Fatalf("failures = %d, want 1", len(result.Failures))
	}
}

func TestColumnEncodesMetadata(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, seriesPath, http.StatusOK, nil)
	sender := NewColumn(newSenderOptions(t, cloud))

	xAxis := "step"
	hidden := true
	system := columnRecord(1, "sys")
	system.GetColumn().ColumnClass = columnv1.ColumnClass_COLUMN_CLASS_SYSTEM
	system.GetColumn().ColumnType = columnv1.ColumnType_COLUMN_TYPE_IMAGE
	system.GetColumn().XAxis = &xAxis
	system.GetColumn().SectionName = "Train"
	system.GetColumn().Hidden = &hidden

	custom := columnRecord(2, "custom")
	custom.GetColumn().ColumnType = columnv1.ColumnType_COLUMN_TYPE_HTML

	result := sender.Upload(context.Background(), []*recordv1.Record{system, custom})
	if len(result.Confirmed) != 2 {
		t.Fatalf("confirmed = %d, want 2", len(result.Confirmed))
	}

	requests := cloud.RequestsFor(http.MethodPost, seriesPath)
	if len(requests) != 1 {
		t.Fatalf("requests = %d, want 1", len(requests))
	}
	series := decodeBody(t, requests[0].Body)["series"].([]any)
	first := series[0].(map[string]any)
	if first["type"] != "SYSTEM" || first["xAxis"] != "step" || first["sectionName"] != "Train" || first["hidden"] != true {
		t.Fatalf("system column = %+v", first)
	}
	second := series[1].(map[string]any)
	if _, ok := second["xAxis"]; ok {
		t.Fatalf("empty xAxis should be omitted: %+v", second)
	}
	if second["type"] != "HTML" {
		t.Fatalf("custom type = %v, want HTML", second["type"])
	}
}

func TestScalarUploadEncodesMetrics(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewScalar(newSenderOptions(t, cloud))

	records := []*recordv1.Record{
		scalarRecord(1, "loss", 3, 0.5),
		scalarRecord(2, "loss", 4, math.NaN()),
		scalarRecord(3, "loss", 5, math.Inf(1)),
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 3 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}

	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(requests) != 1 {
		t.Fatalf("requests = %d, want 1", len(requests))
	}
	body := decodeBody(t, requests[0].Body)
	if body["type"] != "scalar" || body["projectId"] != "proj-1" || body["experimentId"] != "exp-1" {
		t.Fatalf("envelope = %+v", body)
	}
	metrics := body["metrics"].([]any)
	if len(metrics) != 3 {
		t.Fatalf("metrics = %d, want 3", len(metrics))
	}
	first := metrics[0].(map[string]any)
	if first["key"] != "loss" || first["index"] != float64(3) || first["data"] != 0.5 {
		t.Fatalf("first metric = %+v", first)
	}
	if first["create_time"] != "2026-03-10T19:02:03.123400Z" {
		t.Fatalf("create_time = %v", first["create_time"])
	}
	for i := 1; i < 3; i++ {
		if got := metrics[i].(map[string]any)["data"]; got != "nan" {
			t.Fatalf("metric[%d].data = %v, want nan", i, got)
		}
	}
}

func TestScalarMissingTimestampBlocked(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewScalar(newSenderOptions(t, cloud))

	valid := scalarRecord(1, "loss", 1, 1)
	invalid := scalarRecord(2, "loss", 2, 2)
	invalid.GetScalar().Timestamp = nil

	result := sender.Upload(context.Background(), []*recordv1.Record{valid, invalid})
	if len(result.Confirmed) != 2 || len(result.Remaining) != 0 {
		t.Fatalf("result = confirmed %d remaining %d", len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 1 || result.Failures[0].Message != scalarBlockedMessage {
		t.Fatalf("failures = %+v", result.Failures)
	}
	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	metrics := decodeBody(t, requests[0].Body)["metrics"].([]any)
	if len(metrics) != 1 {
		t.Fatalf("metrics = %d, want only the valid record", len(metrics))
	}
}

func TestScalarMissingValueDefaultsToZero(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewScalar(newSenderOptions(t, cloud))

	record := scalarRecord(1, "loss", 1, 1)
	record.GetScalar().Value = nil

	result := sender.Upload(context.Background(), []*recordv1.Record{record})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 0 {
		t.Fatalf("result = %+v, want one confirmed record without failures", result)
	}
	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(requests) != 1 {
		t.Fatalf("metrics requests = %d, want 1", len(requests))
	}
	metrics := decodeBody(t, requests[0].Body)["metrics"].([]any)
	if got := metrics[0].(map[string]any)["data"]; got != float64(0) {
		t.Fatalf("data = %v, want proto default 0", got)
	}
}

func TestScalarRetryableKeepsEncodable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusServiceUnavailable, nil)
	sender := NewScalar(newSenderOptions(t, cloud))

	valid := scalarRecord(1, "loss", 1, 1)
	invalid := scalarRecord(2, "loss", 2, 2)
	invalid.GetScalar().Timestamp = nil

	result := sender.Upload(context.Background(), []*recordv1.Record{valid, invalid})
	if len(result.Confirmed) != 1 || result.Confirmed[0] != invalid {
		t.Fatalf("confirmed = %v, want only the blocked record", result.Confirmed)
	}
	if len(result.Remaining) != 1 || result.Remaining[0] != valid {
		t.Fatalf("remaining = %v, want the encodable record", result.Remaining)
	}
	if len(result.Failures) != 2 {
		t.Fatalf("failures = %d, want transport + blocked", len(result.Failures))
	}
}

func TestScalarDeterministicErrorsTerminate(t *testing.T) {
	cases := []struct {
		name   string
		status int
	}{
		{"house 202", http.StatusAccepted},
		{"unauthorized", http.StatusUnauthorized},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cloud := clienttest.NewCloud(t)
			cloud.ServeJSON(http.MethodPost, metricsPath, tc.status, nil)
			sender := NewScalar(newSenderOptions(t, cloud))

			records := []*recordv1.Record{scalarRecord(1, "loss", 1, 1)}
			result := sender.Upload(context.Background(), records)
			if len(result.Confirmed) != 1 || len(result.Remaining) != 0 || len(result.Failures) != 1 {
				t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
			}
		})
	}
}

func TestContextCancellationKeepsRecords(t *testing.T) {
	t.Run("column", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		sender := NewColumn(newSenderOptions(t, cloud))
		ctx, cancel := context.WithCancel(context.Background())
		cancel()
		result := sender.Upload(ctx, []*recordv1.Record{columnRecord(1, "loss")})
		if len(result.Remaining) != 1 || len(result.Confirmed) != 0 {
			t.Fatalf("result = %+v, want record remaining", result)
		}
		if got := cloud.CountFor(http.MethodPost, seriesPath); got != 0 {
			t.Fatalf("series requests = %d, want 0", got)
		}
	})
	t.Run("scalar", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		sender := NewScalar(newSenderOptions(t, cloud))
		ctx, cancel := context.WithCancel(context.Background())
		cancel()
		result := sender.Upload(ctx, []*recordv1.Record{scalarRecord(1, "loss", 1, 1)})
		if len(result.Remaining) != 1 || len(result.Confirmed) != 0 {
			t.Fatalf("result = %+v, want record remaining", result)
		}
		if got := cloud.CountFor(http.MethodPost, metricsPath); got != 0 {
			t.Fatalf("metrics requests = %d, want 0", got)
		}
	})
	t.Run("log", func(t *testing.T) {
		cloud := clienttest.NewCloud(t)
		sender := NewLog(newSenderOptions(t, cloud))
		ctx, cancel := context.WithCancel(context.Background())
		cancel()
		result := sender.Upload(ctx, []*recordv1.Record{logRecord(1, "line", terminalv1.LogLevel_LOG_LEVEL_INFO, 0)})
		if len(result.Remaining) != 1 || len(result.Confirmed) != 0 {
			t.Fatalf("result = %+v, want record remaining", result)
		}
		if got := cloud.CountFor(http.MethodPost, metricsPath); got != 0 {
			t.Fatalf("metrics requests = %d, want 0", got)
		}
	})
}

func TestLogUploadEncodesMetrics(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewLog(newSenderOptions(t, cloud))

	records := []*recordv1.Record{
		logRecord(1, "hello", terminalv1.LogLevel_LOG_LEVEL_INFO, 7),
		logRecord(2, "oops", terminalv1.LogLevel_LOG_LEVEL_ERROR, 9),
	}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 2 || len(result.Failures) != 0 {
		t.Fatalf("result = confirmed %d failures %d", len(result.Confirmed), len(result.Failures))
	}

	requests := cloud.RequestsFor(http.MethodPost, metricsPath)
	if len(requests) != 1 {
		t.Fatalf("requests = %d, want 1", len(requests))
	}
	body := decodeBody(t, requests[0].Body)
	if body["type"] != "log" {
		t.Fatalf("type = %v, want log", body["type"])
	}
	metrics := body["metrics"].([]any)
	first := metrics[0].(map[string]any)
	if first["level"] != "INFO" || first["message"] != "hello" || first["epoch"] != float64(7) {
		t.Fatalf("first log = %+v", first)
	}
	if first["create_time"] != "2026-03-10T19:02:03.123400Z" {
		t.Fatalf("create_time = %v", first["create_time"])
	}
	if second := metrics[1].(map[string]any); second["level"] != "ERROR" {
		t.Fatalf("second level = %v", second["level"])
	}
}

func TestLogUnknownLevelTerminates(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewLog(newSenderOptions(t, cloud))

	bad := logRecord(1, "weird", terminalv1.LogLevel(7), 0)
	result := sender.Upload(context.Background(), []*recordv1.Record{bad})
	if len(result.Confirmed) != 1 || len(result.Remaining) != 0 {
		t.Fatalf("result = confirmed %d remaining %d", len(result.Confirmed), len(result.Remaining))
	}
	if len(result.Failures) != 1 || result.Failures[0].Message != logUnknownLevelMessage {
		t.Fatalf("failures = %+v", result.Failures)
	}
	if got := cloud.CountFor(http.MethodPost, metricsPath); got != 0 {
		t.Fatalf("metrics requests = %d, want 0", got)
	}
}

func TestLogMissingTimestampBlocked(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusOK, nil)
	sender := NewLog(newSenderOptions(t, cloud))

	bad := logRecord(1, "line", terminalv1.LogLevel_LOG_LEVEL_INFO, 0)
	bad.GetLog().Timestamp = nil
	result := sender.Upload(context.Background(), []*recordv1.Record{bad})
	if len(result.Failures) != 1 || result.Failures[0].Message != logBlockedMessage {
		t.Fatalf("failures = %+v", result.Failures)
	}
}

func TestLogRetryableKeepsEncodable(t *testing.T) {
	cloud := clienttest.NewCloud(t)
	cloud.ServeJSON(http.MethodPost, metricsPath, http.StatusBadGateway, nil)
	sender := NewLog(newSenderOptions(t, cloud))

	records := []*recordv1.Record{logRecord(1, "line", terminalv1.LogLevel_LOG_LEVEL_INFO, 0)}
	result := sender.Upload(context.Background(), records)
	if len(result.Confirmed) != 0 || len(result.Remaining) != 1 || len(result.Failures) != 1 {
		t.Fatalf("result = confirmed %d remaining %d failures %d", len(result.Confirmed), len(result.Remaining), len(result.Failures))
	}
}

func TestScalarBatchMarshalJSONMatchesStdlib(t *testing.T) {
	emptyGot, err := json.Marshal(scalarBatch{})
	if err != nil {
		t.Fatalf("marshal empty scalarBatch: %v", err)
	}
	if string(emptyGot) != "[]" {
		t.Fatalf("empty scalarBatch = %s, want []", emptyGot)
	}

	type refScalarMetric struct {
		Key        string `json:"key"`
		Index      int64  `json:"index"`
		Data       any    `json:"data"`
		CreateTime string `json:"create_time"`
	}

	samples := []struct {
		key   string
		step  int64
		value float64
		secs  int64
		nanos int32
	}{
		{"loss", 0, 0, 1700000000, 0},
		{"neg_zero", 1, math.Copysign(0, -1), 1700000000, 123_000_000},
		{"combine/train_0", 2, 0.5, 1700000001, 123_400_000},
		{"tiny_e5", 3, 1e-5, 1700000002, 123_456_789},
		{"boundary_e6", 4, 1e-6, 1700000003, 5_000_000},
		{"exp_e7", 5, 1e-7, 1700000004, 4_000},
		{"exp_neg_e7", 6, -1e-7, 1700000005, 7},
		{"exp_e10", 7, -1.23456789e-10, 1700000006, 0},
		{"boundary_e20", 8, 1e20, 1700000007, 0},
		{"boundary_e21", 9, 1e21, 1700000008, 0},
		{"max_f64", 10, math.MaxFloat64, 1700000009, 0},
		{"smallest_f64", 11, math.SmallestNonzeroFloat64, 1700000010, 0},
		{`quote"slash\html<>&line` + "\n" + "sep\u2028中文", 12, 42.25, 1700000011, 500_000_000},
		{"nan_val", 13, math.NaN(), 1700000012, 0},
		{"pos_inf", 14, math.Inf(1), 1700000013, 0},
		{"neg_inf", 15, math.Inf(-1), 1700000014, 0},
	}

	batch := make(scalarBatch, 0, len(samples))
	wantList := make([]refScalarMetric, 0, len(samples))
	for _, sample := range samples {
		batch = append(batch, scalarMetric{
			key:   sample.key,
			step:  sample.step,
			value: sample.value,
			secs:  sample.secs,
			nanos: sample.nanos,
		})
		tsStr, ok := encodeTimestamp(&timestamppb.Timestamp{Seconds: sample.secs, Nanos: sample.nanos})
		if !ok {
			t.Fatalf("encodeTimestamp failed for %+v", sample)
		}
		var data any = sample.value
		if math.IsNaN(sample.value) || math.IsInf(sample.value, 0) {
			data = "nan"
		}
		wantList = append(wantList, refScalarMetric{
			Key:        sample.key,
			Index:      sample.step,
			Data:       data,
			CreateTime: tsStr,
		})
	}

	got, err := json.Marshal(batch)
	if err != nil {
		t.Fatalf("marshal scalarBatch: %v", err)
	}
	want, err := json.Marshal(wantList)
	if err != nil {
		t.Fatalf("marshal reference list: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("scalarBatch JSON mismatch:\ngot:  %s\nwant: %s", got, want)
	}
}

func BenchmarkScalarBatchMarshalJSON(b *testing.B) {
	const batchSize = 10_000
	batch := make(scalarBatch, batchSize)
	for i := range batch {
		batch[i] = scalarMetric{
			key:   fmt.Sprintf("combine/train_%d", i%1000),
			step:  int64(i),
			value: float64(i) * 0.123456789,
			secs:  1700000000 + int64(i),
			nanos: 123_456_000,
		}
	}

	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		raw, err := json.Marshal(batch)
		if err != nil {
			b.Fatalf("marshal batch: %v", err)
		}
		if len(raw) == 0 {
			b.Fatal("empty payload")
		}
	}
}
