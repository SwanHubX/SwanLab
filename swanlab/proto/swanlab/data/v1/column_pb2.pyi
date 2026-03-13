from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class ColumnClass(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    COLUMN_CLASS_UNSPECIFIED: _ClassVar[ColumnClass]
    COLUMN_CLASS_CUSTOM: _ClassVar[ColumnClass]

class ColumnType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    COLUMN_TYPE_UNSPECIFIED: _ClassVar[ColumnType]
    COLUMN_TYPE_FLOAT: _ClassVar[ColumnType]
    COLUMN_TYPE_IMAGE: _ClassVar[ColumnType]
    COLUMN_TYPE_AUDIO: _ClassVar[ColumnType]
    COLUMN_TYPE_TEXT: _ClassVar[ColumnType]
    COLUMN_TYPE_ANY: _ClassVar[ColumnType]

class SectionType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    SECTION_TYPE_UNSPECIFIED: _ClassVar[SectionType]
    SECTION_TYPE_PINNED: _ClassVar[SectionType]
    SECTION_TYPE_HIDDEN: _ClassVar[SectionType]
    SECTION_TYPE_PUBLIC: _ClassVar[SectionType]
    SECTION_TYPE_SYSTEM: _ClassVar[SectionType]
COLUMN_CLASS_UNSPECIFIED: ColumnClass
COLUMN_CLASS_CUSTOM: ColumnClass
COLUMN_TYPE_UNSPECIFIED: ColumnType
COLUMN_TYPE_FLOAT: ColumnType
COLUMN_TYPE_IMAGE: ColumnType
COLUMN_TYPE_AUDIO: ColumnType
COLUMN_TYPE_TEXT: ColumnType
COLUMN_TYPE_ANY: ColumnType
SECTION_TYPE_UNSPECIFIED: SectionType
SECTION_TYPE_PINNED: SectionType
SECTION_TYPE_HIDDEN: SectionType
SECTION_TYPE_PUBLIC: SectionType
SECTION_TYPE_SYSTEM: SectionType

class ColumnError(_message.Message):
    __slots__ = ("code", "message")
    CODE_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    code: str
    message: str
    def __init__(self, code: _Optional[str] = ..., message: _Optional[str] = ...) -> None: ...

class YRange(_message.Message):
    __slots__ = ("min", "max")
    MIN_FIELD_NUMBER: _ClassVar[int]
    MAX_FIELD_NUMBER: _ClassVar[int]
    min: float
    max: float
    def __init__(self, min: _Optional[float] = ..., max: _Optional[float] = ...) -> None: ...

class ColumnRecord(_message.Message):
    __slots__ = ("column_class", "column_type", "column_key", "column_name", "column_error", "section_name", "section_type", "y_range", "chart_index", "chart_name", "metric_name", "metric_colors")
    COLUMN_CLASS_FIELD_NUMBER: _ClassVar[int]
    COLUMN_TYPE_FIELD_NUMBER: _ClassVar[int]
    COLUMN_KEY_FIELD_NUMBER: _ClassVar[int]
    COLUMN_NAME_FIELD_NUMBER: _ClassVar[int]
    COLUMN_ERROR_FIELD_NUMBER: _ClassVar[int]
    SECTION_NAME_FIELD_NUMBER: _ClassVar[int]
    SECTION_TYPE_FIELD_NUMBER: _ClassVar[int]
    Y_RANGE_FIELD_NUMBER: _ClassVar[int]
    CHART_INDEX_FIELD_NUMBER: _ClassVar[int]
    CHART_NAME_FIELD_NUMBER: _ClassVar[int]
    METRIC_NAME_FIELD_NUMBER: _ClassVar[int]
    METRIC_COLORS_FIELD_NUMBER: _ClassVar[int]
    column_class: ColumnClass
    column_type: ColumnType
    column_key: str
    column_name: str
    column_error: ColumnError
    section_name: str
    section_type: SectionType
    y_range: YRange
    chart_index: str
    chart_name: str
    metric_name: str
    metric_colors: _containers.RepeatedScalarFieldContainer[str]
    def __init__(self, column_class: _Optional[_Union[ColumnClass, str]] = ..., column_type: _Optional[_Union[ColumnType, str]] = ..., column_key: _Optional[str] = ..., column_name: _Optional[str] = ..., column_error: _Optional[_Union[ColumnError, _Mapping]] = ..., section_name: _Optional[str] = ..., section_type: _Optional[_Union[SectionType, str]] = ..., y_range: _Optional[_Union[YRange, _Mapping]] = ..., chart_index: _Optional[str] = ..., chart_name: _Optional[str] = ..., metric_name: _Optional[str] = ..., metric_colors: _Optional[_Iterable[str]] = ...) -> None: ...
