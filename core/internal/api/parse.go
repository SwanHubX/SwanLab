// @Title        send.go
// @Description  parse the proto message to json
// @Create       cunyue 2025/6/12 18:49

package api //nolint:revive 暂时不改名

import (
	"errors"
	"strings"

	"google.golang.org/protobuf/types/known/structpb"

	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

// Parser provides methods to parse proto messages to JSON sent to the SwanLab Server.
type Parser struct {
}

// ColumnDTO represents the data transfer object for a column record.
// We use RESTFul naming conventions for the fields.
type ColumnDTO struct {
	Class string `json:"class"`
	Type  string `json:"type"`
	Key   string `json:"key"`
	Name  string `json:"name,omitempty"`
	// e.g. {data_class: '', excepted: ''}
	Error       *structpb.Struct `json:"error,omitempty"`
	SectionName string           `json:"section_name,omitempty"`
	SectionType string           `json:"section_type,omitempty"`
	// allow [0, 100] or [0, null] for y_range
	YRange      []*int64 `json:"y_range,omitempty"`
	ChartIndex  string   `json:"chart_index,omitempty"`
	ChartName   string   `json:"chart_name,omitempty"`
	MetricName  string   `json:"metric_name,omitempty"`
	MetricColor []string `json:"metric_color,omitempty"`
}

// ParseColumnRecord parses a ColumnRecord proto message into ColumnDTO.
// Attention: y_range should be serialized to [0, null] in JSON, not ["0", "None"].
func (p *Parser) ParseColumnRecord(record *pb.ColumnRecord) (ColumnDTO, error) {
	// 1. parse the y_range, if error, yRange will be nil
	var yRange []*int64
	if record.GetChartYRange() != nil {
		yRange = []*int64{
			//nolint:protogetter  // We need pointers to represent left null during JSON serialization.
			record.GetChartYRange().Minval,
			//nolint:protogetter  // We need pointers to represent right null during JSON serialization.
			record.GetChartYRange().Maxval,
		}
	}
	// 2. parse the enum type
	// 2.1 column class
	columnClass := strings.TrimPrefix(record.GetColumnClass().String(), "COL_CLASS_")
	// 2.2 column type
	if record.GetColumnType() == pb.ColumnRecord_COL_UNKNOWN {
		return ColumnDTO{}, errors.New("column type is unknown")
	}
	columnType := strings.TrimPrefix(record.GetColumnType().String(), "COL_")
	// 2.3 section type
	sectionType := strings.TrimPrefix(record.GetSectionType().String(), "SEC_")
	// 4. check if the key is empty, if so, return an error
	key := record.GetColumnKey()
	if key == "" {
		return ColumnDTO{}, errors.New("column key cannot be EMPTY")
	}
	// 5. check the metric color length, if it is not empty, it must be 2 elements
	if record.GetMetricColor() != nil && len(record.GetMetricColor()) > 0 && len(record.GetMetricColor()) != 2 {
		return ColumnDTO{}, errors.New("metric color must be empty or have exactly 2 elements")
	}
	return ColumnDTO{
		Class:       columnClass,
		Type:        columnType,
		Key:         record.GetColumnKey(),
		Name:        record.GetColumnName(),
		Error:       record.GetColumnError(),
		SectionName: record.GetSectionName(),
		SectionType: sectionType,
		YRange:      yRange,
		ChartIndex:  record.GetChartIndex(),
		ChartName:   record.GetChartName(),
		MetricName:  record.GetMetricName(),
		MetricColor: record.GetMetricColor(),
	}, nil
}
