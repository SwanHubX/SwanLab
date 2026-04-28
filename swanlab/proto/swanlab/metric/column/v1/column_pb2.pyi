from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class ColumnType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    COLUMN_TYPE_UNSPECIFIED: _ClassVar[ColumnType]
    COLUMN_TYPE_SCALAR: _ClassVar[ColumnType]
    COLUMN_TYPE_IMAGE: _ClassVar[ColumnType]
    COLUMN_TYPE_AUDIO: _ClassVar[ColumnType]
    COLUMN_TYPE_TEXT: _ClassVar[ColumnType]
    COLUMN_TYPE_VIDEO: _ClassVar[ColumnType]
    COLUMN_TYPE_ECHARTS: _ClassVar[ColumnType]
    COLUMN_TYPE_OBJECT3D: _ClassVar[ColumnType]
    COLUMN_TYPE_MOLECULE: _ClassVar[ColumnType]

class ColumnClass(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    COLUMN_CLASS_UNSPECIFIED: _ClassVar[ColumnClass]
    COLUMN_CLASS_CUSTOM: _ClassVar[ColumnClass]
    COLUMN_CLASS_SYSTEM: _ClassVar[ColumnClass]

class ChartType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    CHART_TYPE_UNSPECIFIED: _ClassVar[ChartType]
    CHART_TYPE_LINE: _ClassVar[ChartType]
    CHART_TYPE_BAR: _ClassVar[ChartType]
    CHART_TYPE_SCALAR: _ClassVar[ChartType]
    CHART_TYPE_IMAGE: _ClassVar[ChartType]
    CHART_TYPE_AUDIO: _ClassVar[ChartType]
    CHART_TYPE_TEXT: _ClassVar[ChartType]
    CHART_TYPE_VIDEO: _ClassVar[ChartType]
    CHART_TYPE_ECHARTS: _ClassVar[ChartType]
    CHART_TYPE_OBJECT3D: _ClassVar[ChartType]
    CHART_TYPE_MOLECULE: _ClassVar[ChartType]

class SectionType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    SECTION_TYPE_UNSPECIFIED: _ClassVar[SectionType]
    SECTION_TYPE_PINNED: _ClassVar[SectionType]
    SECTION_TYPE_HIDDEN: _ClassVar[SectionType]
    SECTION_TYPE_PUBLIC: _ClassVar[SectionType]
    SECTION_TYPE_SYSTEM: _ClassVar[SectionType]
COLUMN_TYPE_UNSPECIFIED: ColumnType
COLUMN_TYPE_SCALAR: ColumnType
COLUMN_TYPE_IMAGE: ColumnType
COLUMN_TYPE_AUDIO: ColumnType
COLUMN_TYPE_TEXT: ColumnType
COLUMN_TYPE_VIDEO: ColumnType
COLUMN_TYPE_ECHARTS: ColumnType
COLUMN_TYPE_OBJECT3D: ColumnType
COLUMN_TYPE_MOLECULE: ColumnType
COLUMN_CLASS_UNSPECIFIED: ColumnClass
COLUMN_CLASS_CUSTOM: ColumnClass
COLUMN_CLASS_SYSTEM: ColumnClass
CHART_TYPE_UNSPECIFIED: ChartType
CHART_TYPE_LINE: ChartType
CHART_TYPE_BAR: ChartType
CHART_TYPE_SCALAR: ChartType
CHART_TYPE_IMAGE: ChartType
CHART_TYPE_AUDIO: ChartType
CHART_TYPE_TEXT: ChartType
CHART_TYPE_VIDEO: ChartType
CHART_TYPE_ECHARTS: ChartType
CHART_TYPE_OBJECT3D: ChartType
CHART_TYPE_MOLECULE: ChartType
SECTION_TYPE_UNSPECIFIED: SectionType
SECTION_TYPE_PINNED: SectionType
SECTION_TYPE_HIDDEN: SectionType
SECTION_TYPE_PUBLIC: SectionType
SECTION_TYPE_SYSTEM: SectionType

class YRange(_message.Message):
    __slots__ = ("min", "max")
    MIN_FIELD_NUMBER: _ClassVar[int]
    MAX_FIELD_NUMBER: _ClassVar[int]
    min: float
    max: float
    def __init__(self, min: _Optional[float] = ..., max: _Optional[float] = ...) -> None: ...

class MetricColors(_message.Message):
    __slots__ = ("light", "dark")
    LIGHT_FIELD_NUMBER: _ClassVar[int]
    DARK_FIELD_NUMBER: _ClassVar[int]
    light: str
    dark: str
    def __init__(self, light: _Optional[str] = ..., dark: _Optional[str] = ...) -> None: ...

class ColumnRecord(_message.Message):
    __slots__ = ("column_class", "column_type", "column_key", "column_name", "section_name", "section_type", "y_range", "chart_index", "chart_name", "chart_type", "metric_name", "metric_colors")
    COLUMN_CLASS_FIELD_NUMBER: _ClassVar[int]
    COLUMN_TYPE_FIELD_NUMBER: _ClassVar[int]
    COLUMN_KEY_FIELD_NUMBER: _ClassVar[int]
    COLUMN_NAME_FIELD_NUMBER: _ClassVar[int]
    SECTION_NAME_FIELD_NUMBER: _ClassVar[int]
    SECTION_TYPE_FIELD_NUMBER: _ClassVar[int]
    Y_RANGE_FIELD_NUMBER: _ClassVar[int]
    CHART_INDEX_FIELD_NUMBER: _ClassVar[int]
    CHART_NAME_FIELD_NUMBER: _ClassVar[int]
    CHART_TYPE_FIELD_NUMBER: _ClassVar[int]
    METRIC_NAME_FIELD_NUMBER: _ClassVar[int]
    METRIC_COLORS_FIELD_NUMBER: _ClassVar[int]
    column_class: ColumnClass
    column_type: ColumnType
    column_key: str
    column_name: str
    section_name: str
    section_type: SectionType
    y_range: YRange
    chart_index: str
    chart_name: str
    chart_type: ChartType
    metric_name: str
    metric_colors: MetricColors
    def __init__(self, column_class: _Optional[_Union[ColumnClass, str]] = ..., column_type: _Optional[_Union[ColumnType, str]] = ..., column_key: _Optional[str] = ..., column_name: _Optional[str] = ..., section_name: _Optional[str] = ..., section_type: _Optional[_Union[SectionType, str]] = ..., y_range: _Optional[_Union[YRange, _Mapping]] = ..., chart_index: _Optional[str] = ..., chart_name: _Optional[str] = ..., chart_type: _Optional[_Union[ChartType, str]] = ..., metric_name: _Optional[str] = ..., metric_colors: _Optional[_Union[MetricColors, _Mapping]] = ...) -> None: ...
