// @Title        send.go
// @Description  parse the proto message to json
// @Create       cunyue 2025/6/12 18:49

package api

import (
	"errors"
	"strconv"

	"github.com/SwanHubX/SwanLab/core/pkg/pb"
)

// Some constants for the parser.
const (
	minRangeLen = 2
)

// Parser provides methods to parse proto messages to JSON sent to the SwanLab Server.
type Parser struct {
}

// convertRange will convert the YRange from a proto message (like ["0", "None"]) to a JSON-compatible format(like [0, null]).
// We expect the input to be a slice of strings with 2 elements, any other length will return an empty slice.
func (p *Parser) convertRange(r []string) []interface{} {
	var left, right interface{} = nil, nil
	// must be at least 2 elements
	if len(r) != minRangeLen {
		return []interface{}{}
	}
	if val, err := strconv.Atoi(r[0]); err == nil {
		left = val
	}
	if val, err := strconv.Atoi(r[1]); err == nil {
		right = val
	}
	return []interface{}{left, right}
}

// validateColumnClass will validate the columnClass.
func (p *Parser) validateColumnClass(record *pb.ColumnRecord) (string, error) {
	var class string
	switch record.GetColumnClass() {
	case pb.ColumnRecord_COL_CUSTOM:
		class = "CUSTOM"
	case pb.ColumnRecord_COL_SYSTEM:
		class = "SYSTEM"
	default:
		return "", errors.New("invalid column class: " + record.GetColumnClass().String())
	}
	return class, nil
}

// validateColumnType will validate the columnType.
func (p *Parser) validateColumnType(record *pb.ColumnRecord) (string, error) {
	switch record.GetColumnType() {
	case pb.ColumnRecord_COL_FLOAT:
		return "FLOAT", nil
	case pb.ColumnRecord_COL_IMAGE:
		return "IMAGE", nil
	case pb.ColumnRecord_COL_AUDIO:
		return "AUDIO", nil
	case pb.ColumnRecord_COL_TEXT:
		return "TEXT", nil
	case pb.ColumnRecord_COL_OBJECT3D:
		return "OBJECT3D", nil
	case pb.ColumnRecord_COL_MOLECULE:
		return "MOLECULE", nil
	case pb.ColumnRecord_COL_ECHARTS:
		return "ECHARTS", nil
	default:
		return "", errors.New("invalid column type: " + record.GetColumnType().String())
	}
}

// validateSectionType will validate the sectionType.
func (p *Parser) validateSectionType(record *pb.ColumnRecord) (string, error) {
	var sectionType string
	switch record.GetSectionType() {
	case pb.ColumnRecord_SEC_CUSTOM:
		sectionType = "CUSTOM"
	case pb.ColumnRecord_SEC_SYSTEM:
		sectionType = "SYSTEM"
	case pb.ColumnRecord_SEC_PINNED:
		sectionType = "PINNED"
	case pb.ColumnRecord_SEC_HIDDEN:
		sectionType = "HIDDEN"
	case pb.ColumnRecord_SEC_PUBLIC:
		sectionType = "PUBLIC"
	default:
		return "", errors.New("invalid section type: " + record.GetSectionType().String())
	}
	return sectionType, nil
}

type ColumnDTO struct {
	Class       string      `json:"class"`
	Type        string      `json:"type"`
	Key         string      `json:"key"`
	Name        string      `json:"name,omitempty"`
	Error       interface{} `json:"error,omitempty"`
	SectionName string      `json:"section_name,omitempty"`
	SectionType string      `json:"section_type,omitempty"`
	// allow [0, 100] or [0, null] for y_range
	YRange      []interface{} `json:"y_range,omitempty"`
	ChartIndex  string        `json:"chart_index,omitempty"`
	ChartName   string        `json:"chart_name,omitempty"`
	MetricName  string        `json:"metric_name,omitempty"`
	MetricColor []string      `json:"metric_color,omitempty"`
}

// ParseColumnRecord parses a ColumnRecord proto message into a map[string]interface{} for JSON serialization.
// Attention: y_range should be serialized to [0, null] in JSON, not ["0", "None"].
func (p *Parser) ParseColumnRecord(record *pb.ColumnRecord) (ColumnDTO, error) {
	// 1. parse the y_range, if error, yRange will be nil
	yRange := p.convertRange(record.GetChartYRange())
	// 2. parse the enum type
	// 2.1 column class
	class, err := p.validateColumnClass(record)
	if err != nil {
		return ColumnDTO{}, err
	}
	// 2.2 column type
	columnType, err := p.validateColumnType(record)
	if err != nil {
		return ColumnDTO{}, err
	}
	// 2.3 section type
	sectionType, err := p.validateSectionType(record)
	if err != nil {
		return ColumnDTO{}, err
	}
	// 4. check if the key is empty, if so, return an error
	key := record.GetColumnKey()
	if key == "" {
		return ColumnDTO{}, errors.New("column key cannot be EMPTY")
	}
	// 5. check the metric color length, if it is not empty, it must be 2 elements
	if len(record.GetMetricColor()) > 0 && len(record.GetMetricColor()) != 2 {
		return ColumnDTO{}, errors.New("metric color must be empty or have exactly 2 elements")
	}
	return ColumnDTO{
		Class:       class,
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
