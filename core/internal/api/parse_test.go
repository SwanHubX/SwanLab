// @Title        parse_test.go
// @Description  test the parse function
// @Create       cunyue 2025/6/12 18:50

package api_test

import (
	"encoding/json"
	"math/rand"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/SwanHubX/SwanLab/core/internal/api"
	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

// TestRangeSerialization tests the ParseColumnRecord function of the Parser.
// Just a simple test to ensure that the yRange can be parsed correctly.
func TestRangeSerialization(t *testing.T) {
	// Define a test case with a valid yRange
	tests := []struct {
		name      string
		input     []string
		expected  []interface{}
		yRangeStr string
	}{
		{
			name:      "Valid yRange",
			input:     []string{"1", "2"},
			expected:  []interface{}{1, 2},
			yRangeStr: "[1,2]",
		},
		{
			name:      "Too many values in yRange",
			input:     []string{"1", "2", "3"},
			expected:  []interface{}{},
			yRangeStr: "[]",
		},
		{
			name:      "Left is none",
			input:     []string{"None", "2"},
			expected:  []interface{}{nil, 2},
			yRangeStr: "[null,2]",
		},
		{
			name:      "Right is none",
			input:     []string{"1", "None"},
			expected:  []interface{}{1, nil},
			yRangeStr: "[1,null]",
		},
		{
			name:     "Invalid yRange",
			input:    []string{"x", "invalid"},
			expected: []interface{}{nil, nil},
		},
	}
	parser := api.Parser{}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			record := &pb.ColumnRecord{
				ColumnKey:   "test-key",
				ColumnClass: pb.ColumnRecord_COL_CUSTOM,
				ColumnType:  pb.ColumnRecord_COL_ECHARTS,
				SectionName: "",
				SectionType: pb.ColumnRecord_SEC_PUBLIC,
				ChartName:   "test-chart-name",
				ChartYRange: tt.input,
				ChartIndex:  "test-chart-index" + string(rune(rand.Intn(1000))),
				MetricName:  "test-metric-name",
				MetricColor: []string{"#FF0000", "#0000"},
			}
			// Call the ParseColumnRecord function with the test input
			result, err := parser.ParseColumnRecord(record)
			if err != nil {
				t.Errorf("Parse error: %s", err)
				return
			}
			// Check base fields
			assert.Equal(t, record.GetColumnKey(), result.Key)
			assert.Empty(t, result.Name)
			assert.Equal(t, "CUSTOM", result.Class)
			assert.Equal(t, "ECHARTS", result.Type)
			assert.Empty(t, result.SectionName)
			assert.Equal(t, "PUBLIC", result.SectionType)
			assert.Equal(t, record.GetChartName(), result.ChartName)
			assert.Equal(t, record.GetChartIndex(), result.ChartIndex)
			assert.Equal(t, record.GetMetricName(), result.MetricName)
			assert.Equal(t, record.GetMetricColor()[0], result.MetricColor[0])
			assert.Equal(t, record.GetMetricColor()[1], result.MetricColor[1])
			// Check if the yRange is parsed correctly
			if len(result.YRange) != len(tt.expected) {
				t.Errorf("Expected yRange length %d, got %d", len(tt.expected), len(result.YRange))
			}
			// check if the yRange values are equal
			if len(tt.expected) == 2 {
				assert.Equal(t, tt.expected[0], result.YRange[0], "Left value mismatch")
				assert.Equal(t, tt.expected[1], result.YRange[1], "Right value mismatch")
			} else {
				assert.Empty(t, result.YRange, "Expected empty yRange for invalid input")
			}
			// JSON serialization, check if it can be serialized without error
			// Also check the yRange is serialized correctly
			jsonData, err := json.Marshal(result)
			require.NoError(t, err)
			jsonStr := string(jsonData)
			if len(tt.expected) == 2 {
				require.Contains(t, jsonStr, tt.yRangeStr, "Expected yRange to be serialized correctly")
			} else {
				require.NotContains(t, "y_range", jsonStr, "Expected yRange to be empty for invalid input")
			}
		})
	}
}
