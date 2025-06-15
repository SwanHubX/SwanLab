// @Title        parse_test.go
// @Description  test the parse function
// @Create       cunyue 2025/6/12 18:50

package api_test

import (
	"encoding/json"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"

	"github.com/SwanHubX/SwanLab/core/internal/api"
	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

// TestParser_ParseColumnRecord tests the basic function of the ParseColumnRecord function in the Parser.
func TestParser_ParseColumnRecord(t *testing.T) {
	record := pb.ColumnRecord{
		ColumnKey:   "test-key",
		SectionName: "",
		ChartName:   "test-chart-name",
		ChartYRange: nil,
		ColumnType:  pb.ColumnRecord_COL_FLOAT,
		ChartIndex:  "test-chart-index",
		MetricName:  "test-metric-name",
		MetricColor: []string{"#FF0000", "#0000"},
	}
	parser := api.Parser{}
	columnDTO, err := parser.ParseColumnRecord(&record)
	require.NoError(t, err)
	assert.Equal(t, "test-key", columnDTO.Key)
	assert.Equal(t, "CUSTOM", columnDTO.Class)
	assert.Equal(t, "FLOAT", columnDTO.Type)
	assert.Equal(t, "PUBLIC", columnDTO.SectionType)
	assert.Empty(t, columnDTO.SectionName)
	assert.Equal(t, "test-chart-name", columnDTO.ChartName)
	assert.Equal(t, "test-chart-index", columnDTO.ChartIndex)
	assert.Equal(t, "test-metric-name", columnDTO.MetricName)
	assert.Nil(t, columnDTO.YRange, "Expected YRange to be nil when ChartYRange is not set")
	assert.Equal(t, "#FF0000", columnDTO.MetricColor[0], "Expected first metric color to be #FF0000")
	assert.Equal(t, "#0000", columnDTO.MetricColor[1], "Expected second metric color to be #0000")
}

// TestParser_ParseColumnRecord_YRange tests the ParseColumnRecord function of the Parser.
// This test is specifically for the yRange parsing logic, and the json serialization of the yRange field.
func TestParser_ParseColumnRecord_YRange(t *testing.T) {
	minval := int64(1)
	maxval := int64(2)
	tests := []struct {
		name     string
		input    *pb.Range
		expected *[]*int64
		// '' indicates that yRange should not be in the JSON string
		jsonStr string
	}{
		{
			name:     "Valid yRange",
			input:    &pb.Range{Minval: &minval, Maxval: &maxval},
			expected: &[]*int64{&minval, &maxval},
			jsonStr:  "[1,2]",
		},
		{
			name:     "Left is None",
			input:    &pb.Range{Minval: nil, Maxval: &maxval},
			expected: &[]*int64{nil, &maxval},
			jsonStr:  "[null,2]",
		},
		{
			name:     "Right is None",
			input:    &pb.Range{Minval: &minval, Maxval: nil},
			expected: &[]*int64{&minval, nil},
			jsonStr:  "[1,null]",
		},
		{
			name:     "Both are None",
			input:    &pb.Range{Minval: nil, Maxval: nil},
			expected: &[]*int64{nil, nil},
			jsonStr:  "[null,null]",
		},
		{
			name:     "Nil yRange",
			input:    nil,
			expected: nil,
		},
	}
	parser := api.Parser{}
	for _, tt := range tests {
		t.Run(
			tt.name, func(t *testing.T) {
				record := &pb.ColumnRecord{
					ColumnKey:   "test-y-range",
					ColumnType:  pb.ColumnRecord_COL_FLOAT,
					ChartYRange: tt.input,
				}
				columnDto, err := parser.ParseColumnRecord(record)
				require.NoError(t, err, "ParseColumnRecord should not return an error")
				if tt.expected != nil {
					// yRange has been set
					assert.Equal(t, (*tt.expected)[0], columnDto.YRange[0], "Left value mismatch")
					assert.Equal(t, (*tt.expected)[1], columnDto.YRange[1], "Right value mismatch")
				} else {
					// yRange is nil
					assert.Nil(t, columnDto.YRange, "Expected YRange to be nil")
				}
				// Test JSON serialization
				jsonData, err := json.Marshal(columnDto)
				require.NoError(t, err, "JSON serialization should not return an error")
				jsonStr := string(jsonData)
				if tt.jsonStr != "" {
					require.Contains(t, jsonStr, tt.jsonStr, "Expected yRange to be serialized correctly")
				} else {
					require.NotContains(t, jsonStr, "y_range", "Expected yRange to be empty for invalid input")
				}
			},
		)
	}
}

func TestParser_ParseColumnRecord_ColumnClass(t *testing.T) {
	tests := []struct {
		name     string
		input    pb.ColumnRecord_ColumnClass
		expected string
	}{
		{
			name:     "Custom Column Class",
			input:    pb.ColumnRecord_COL_CLASS_CUSTOM,
			expected: "CUSTOM",
		},
		{
			name:     "System Column Class",
			input:    pb.ColumnRecord_COL_CLASS_SYSTEM,
			expected: "SYSTEM",
		},
		{
			name:     "Default Column Class",
			expected: "CUSTOM",
		},
	}

	parser := api.Parser{}
	for _, tt := range tests {
		t.Run(
			tt.name, func(t *testing.T) {
				record := &pb.ColumnRecord{
					ColumnKey:   "test-key",
					ColumnType:  pb.ColumnRecord_COL_FLOAT,
					ColumnClass: tt.input,
				}
				columnDto, err := parser.ParseColumnRecord(record)
				require.NoError(t, err, "ParseColumnRecord should not return an error")
				assert.Equal(t, tt.expected, columnDto.Class, "Column class mismatch")
			},
		)
	}
}

func TestParser_ParseColumnRecord_ColumnType(t *testing.T) {
	tests := []struct {
		name     string
		input    pb.ColumnRecord_ColumnType
		expected string
		err      bool
	}{
		{
			name:     "FLOAT Column Type",
			input:    pb.ColumnRecord_COL_FLOAT,
			expected: "FLOAT",
		},
		{
			name:     "Image Column Type",
			input:    pb.ColumnRecord_COL_IMAGE,
			expected: "IMAGE",
		},
		{
			name:     "Audio Column Type",
			input:    pb.ColumnRecord_COL_AUDIO,
			expected: "AUDIO",
		},
		{
			name:     "Text Column Type",
			input:    pb.ColumnRecord_COL_TEXT,
			expected: "TEXT",
		},
		{
			name:     "Object3D Column Type",
			input:    pb.ColumnRecord_COL_OBJECT3D,
			expected: "OBJECT3D",
		},
		{
			name:     "Molecule Column Type",
			input:    pb.ColumnRecord_COL_MOLECULE,
			expected: "MOLECULE",
		},
		{
			name:     "ECharts Column Type",
			input:    pb.ColumnRecord_COL_ECHARTS,
			expected: "ECHARTS",
		},
		{
			name: "Unknown Column Type",
			err:  true,
		},
	}
	parser := api.Parser{}
	for _, tt := range tests {
		t.Run(
			tt.name, func(t *testing.T) {
				record := &pb.ColumnRecord{
					ColumnKey:  "test-key",
					ColumnType: tt.input,
				}
				columnDto, err := parser.ParseColumnRecord(record)
				if tt.err {
					assert.Error(t, err, "Expected an error for invalid column type")
					return
				}
				require.NoError(t, err, "ParseColumnRecord should not return an error")
				assert.Equal(t, tt.expected, columnDto.Type, "Column type mismatch")
			},
		)
	}
}

func TestParser_ParseColumnRecord_SectionType(t *testing.T) {
	tests := []struct {
		name     string
		input    pb.ColumnRecord_SectionType
		expected string
	}{
		{
			name:     "Public Section Type",
			input:    pb.ColumnRecord_SEC_PUBLIC,
			expected: "PUBLIC",
		},
		{
			name:     "System Section Type",
			input:    pb.ColumnRecord_SEC_SYSTEM,
			expected: "SYSTEM",
		},
		{
			name:     "Default Section Type",
			expected: "PUBLIC", // Default value if not set
		},
		{
			name:     "Pinned Section Type",
			input:    pb.ColumnRecord_SEC_PINNED,
			expected: "PINNED",
		},
		{
			name:     "Hidden Section Type",
			input:    pb.ColumnRecord_SEC_HIDDEN,
			expected: "HIDDEN",
		},
	}
	parser := api.Parser{}
	for _, tt := range tests {
		t.Run(
			tt.name, func(t *testing.T) {
				record := &pb.ColumnRecord{
					ColumnKey:   "test-key",
					ColumnType:  pb.ColumnRecord_COL_FLOAT,
					SectionType: tt.input,
				}
				columnDto, err := parser.ParseColumnRecord(record)
				require.NoError(t, err, "ParseColumnRecord should not return an error")
				assert.Equal(t, tt.expected, columnDto.SectionType, "Section type mismatch")
			},
		)
	}
}

func TestParser_ParseColumnRecord_Key(t *testing.T) {
	tests := []struct {
		name     string
		input    string
		expected string
		err      bool
	}{
		{
			name:     "Valid Key",
			input:    "valid-key",
			expected: "valid-key",
		},
		{
			name: "Empty Key",
			err:  true,
		},
	}

	parser := api.Parser{}
	for _, tt := range tests {
		t.Run(
			tt.name, func(t *testing.T) {
				record := &pb.ColumnRecord{
					ColumnKey:  tt.input,
					ColumnType: pb.ColumnRecord_COL_FLOAT,
				}
				columnDto, err := parser.ParseColumnRecord(record)
				if tt.err {
					assert.Error(t, err, "Expected an error for empty key")
					return
				}
				require.NoError(t, err, "ParseColumnRecord should not return an error")
				assert.Equal(t, tt.expected, columnDto.Key, "Column key mismatch")
			},
		)
	}
}

func TestParser_ParseColumnRecord_MetricColor(t *testing.T) {
	tests := []struct {
		name     string
		input    []string
		expected []string
		err      bool
	}{
		{
			name:     "Valid Metric Color",
			input:    []string{"#FF0000", "#00FF00"},
			expected: []string{"#FF0000", "#00FF00"},
		},
		{
			name:  "Empty Metric Color",
			input: nil,
			err:   false, // Should not return an error for empty metric color
		},
		{
			name:  "Invalid Metric Color Length",
			input: []string{"#FF0000", "#00FF00", "#0000FF"},
			err:   true, // Should return an error for invalid length
		},
		{
			name:  "Single Metric Color",
			input: []string{"#FF0000"},
			err:   true, // Should return an error for single color
		},
	}

	parser := api.Parser{}
	for _, tt := range tests {
		t.Run(
			tt.name, func(t *testing.T) {
				record := &pb.ColumnRecord{
					ColumnKey:   "test-metric-color",
					ColumnType:  pb.ColumnRecord_COL_FLOAT,
					MetricColor: tt.input,
				}
				columnDto, err := parser.ParseColumnRecord(record)
				if tt.err {
					assert.Error(t, err, "Expected an error for invalid metric color length")
					return
				}
				require.NoError(t, err, "ParseColumnRecord should not return an error")
				assert.Equal(t, tt.expected, columnDto.MetricColor, "Metric color mismatch")
			},
		)
	}
}
