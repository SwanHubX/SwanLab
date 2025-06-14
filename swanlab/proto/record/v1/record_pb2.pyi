from google.protobuf import struct_pb2 as _struct_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class SetupRecord(_message.Message):
    __slots__ = ("name", "workspace", "public", "experiment_name", "experiment_description", "experiment_tags", "start_time")
    NAME_FIELD_NUMBER: _ClassVar[int]
    WORKSPACE_FIELD_NUMBER: _ClassVar[int]
    PUBLIC_FIELD_NUMBER: _ClassVar[int]
    EXPERIMENT_NAME_FIELD_NUMBER: _ClassVar[int]
    EXPERIMENT_DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    EXPERIMENT_TAGS_FIELD_NUMBER: _ClassVar[int]
    START_TIME_FIELD_NUMBER: _ClassVar[int]
    name: str
    workspace: str
    public: bool
    experiment_name: str
    experiment_description: str
    experiment_tags: _containers.RepeatedScalarFieldContainer[str]
    start_time: str
    def __init__(self, name: _Optional[str] = ..., workspace: _Optional[str] = ..., public: bool = ..., experiment_name: _Optional[str] = ..., experiment_description: _Optional[str] = ..., experiment_tags: _Optional[_Iterable[str]] = ..., start_time: _Optional[str] = ...) -> None: ...

class TeardownRecord(_message.Message):
    __slots__ = ("state", "error_message", "end_time")
    class State(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        STATE_SUCCESS: _ClassVar[TeardownRecord.State]
        STATE_FAILED: _ClassVar[TeardownRecord.State]
    STATE_SUCCESS: TeardownRecord.State
    STATE_FAILED: TeardownRecord.State
    STATE_FIELD_NUMBER: _ClassVar[int]
    ERROR_MESSAGE_FIELD_NUMBER: _ClassVar[int]
    END_TIME_FIELD_NUMBER: _ClassVar[int]
    state: TeardownRecord.State
    error_message: str
    end_time: str
    def __init__(self, state: _Optional[_Union[TeardownRecord.State, str]] = ..., error_message: _Optional[str] = ..., end_time: _Optional[str] = ...) -> None: ...

class RuntimeRecord(_message.Message):
    __slots__ = ("conda_filename", "pip_filename", "config_filename", "metadata_filename")
    CONDA_FILENAME_FIELD_NUMBER: _ClassVar[int]
    PIP_FILENAME_FIELD_NUMBER: _ClassVar[int]
    CONFIG_FILENAME_FIELD_NUMBER: _ClassVar[int]
    METADATA_FILENAME_FIELD_NUMBER: _ClassVar[int]
    conda_filename: str
    pip_filename: str
    config_filename: str
    metadata_filename: str
    def __init__(self, conda_filename: _Optional[str] = ..., pip_filename: _Optional[str] = ..., config_filename: _Optional[str] = ..., metadata_filename: _Optional[str] = ...) -> None: ...

class Range(_message.Message):
    __slots__ = ("minval", "maxval")
    MINVAL_FIELD_NUMBER: _ClassVar[int]
    MAXVAL_FIELD_NUMBER: _ClassVar[int]
    minval: int
    maxval: int
    def __init__(self, minval: _Optional[int] = ..., maxval: _Optional[int] = ...) -> None: ...

class ColumnRecord(_message.Message):
    __slots__ = ("column_key", "column_name", "column_class", "column_type", "column_error", "section_name", "section_type", "chart_name", "chart_index", "chart_y_range", "metric_name", "metric_color")
    class ColumClass(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        COL_CUSTOM: _ClassVar[ColumnRecord.ColumClass]
        COL_SYSTEM: _ClassVar[ColumnRecord.ColumClass]
    COL_CUSTOM: ColumnRecord.ColumClass
    COL_SYSTEM: ColumnRecord.ColumClass
    class ColumnType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        COL_FLOAT: _ClassVar[ColumnRecord.ColumnType]
        COL_IMAGE: _ClassVar[ColumnRecord.ColumnType]
        COL_AUDIO: _ClassVar[ColumnRecord.ColumnType]
        COL_TEXT: _ClassVar[ColumnRecord.ColumnType]
        COL_OBJECT3D: _ClassVar[ColumnRecord.ColumnType]
        COL_MOLECULE: _ClassVar[ColumnRecord.ColumnType]
        COL_ECHARTS: _ClassVar[ColumnRecord.ColumnType]
    COL_FLOAT: ColumnRecord.ColumnType
    COL_IMAGE: ColumnRecord.ColumnType
    COL_AUDIO: ColumnRecord.ColumnType
    COL_TEXT: ColumnRecord.ColumnType
    COL_OBJECT3D: ColumnRecord.ColumnType
    COL_MOLECULE: ColumnRecord.ColumnType
    COL_ECHARTS: ColumnRecord.ColumnType
    class SectionType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        SEC_PINNED: _ClassVar[ColumnRecord.SectionType]
        SEC_HIDDEN: _ClassVar[ColumnRecord.SectionType]
        SEC_PUBLIC: _ClassVar[ColumnRecord.SectionType]
        SEC_SYSTEM: _ClassVar[ColumnRecord.SectionType]
        SEC_CUSTOM: _ClassVar[ColumnRecord.SectionType]
    SEC_PINNED: ColumnRecord.SectionType
    SEC_HIDDEN: ColumnRecord.SectionType
    SEC_PUBLIC: ColumnRecord.SectionType
    SEC_SYSTEM: ColumnRecord.SectionType
    SEC_CUSTOM: ColumnRecord.SectionType
    COLUMN_KEY_FIELD_NUMBER: _ClassVar[int]
    COLUMN_NAME_FIELD_NUMBER: _ClassVar[int]
    COLUMN_CLASS_FIELD_NUMBER: _ClassVar[int]
    COLUMN_TYPE_FIELD_NUMBER: _ClassVar[int]
    COLUMN_ERROR_FIELD_NUMBER: _ClassVar[int]
    SECTION_NAME_FIELD_NUMBER: _ClassVar[int]
    SECTION_TYPE_FIELD_NUMBER: _ClassVar[int]
    CHART_NAME_FIELD_NUMBER: _ClassVar[int]
    CHART_INDEX_FIELD_NUMBER: _ClassVar[int]
    CHART_Y_RANGE_FIELD_NUMBER: _ClassVar[int]
    METRIC_NAME_FIELD_NUMBER: _ClassVar[int]
    METRIC_COLOR_FIELD_NUMBER: _ClassVar[int]
    column_key: str
    column_name: str
    column_class: ColumnRecord.ColumClass
    column_type: ColumnRecord.ColumnType
    column_error: _struct_pb2.Struct
    section_name: str
    section_type: ColumnRecord.SectionType
    chart_name: str
    chart_index: str
    chart_y_range: Range
    metric_name: str
    metric_color: _containers.RepeatedScalarFieldContainer[str]
    def __init__(self, column_key: _Optional[str] = ..., column_name: _Optional[str] = ..., column_class: _Optional[_Union[ColumnRecord.ColumClass, str]] = ..., column_type: _Optional[_Union[ColumnRecord.ColumnType, str]] = ..., column_error: _Optional[_Union[_struct_pb2.Struct, _Mapping]] = ..., section_name: _Optional[str] = ..., section_type: _Optional[_Union[ColumnRecord.SectionType, str]] = ..., chart_name: _Optional[str] = ..., chart_index: _Optional[str] = ..., chart_y_range: _Optional[_Union[Range, _Mapping]] = ..., metric_name: _Optional[str] = ..., metric_color: _Optional[_Iterable[str]] = ...) -> None: ...

class MediaRecord(_message.Message):
    __slots__ = ("index", "epoch", "create_time", "key", "key_encoded", "kid", "data", "more")
    INDEX_FIELD_NUMBER: _ClassVar[int]
    EPOCH_FIELD_NUMBER: _ClassVar[int]
    CREATE_TIME_FIELD_NUMBER: _ClassVar[int]
    KEY_FIELD_NUMBER: _ClassVar[int]
    KEY_ENCODED_FIELD_NUMBER: _ClassVar[int]
    KID_FIELD_NUMBER: _ClassVar[int]
    DATA_FIELD_NUMBER: _ClassVar[int]
    MORE_FIELD_NUMBER: _ClassVar[int]
    index: str
    epoch: str
    create_time: str
    key: str
    key_encoded: str
    kid: str
    data: _containers.RepeatedScalarFieldContainer[str]
    more: _containers.RepeatedScalarFieldContainer[str]
    def __init__(self, index: _Optional[str] = ..., epoch: _Optional[str] = ..., create_time: _Optional[str] = ..., key: _Optional[str] = ..., key_encoded: _Optional[str] = ..., kid: _Optional[str] = ..., data: _Optional[_Iterable[str]] = ..., more: _Optional[_Iterable[str]] = ...) -> None: ...

class ScalarRecord(_message.Message):
    __slots__ = ("index", "epoch", "create_time", "key", "data")
    INDEX_FIELD_NUMBER: _ClassVar[int]
    EPOCH_FIELD_NUMBER: _ClassVar[int]
    CREATE_TIME_FIELD_NUMBER: _ClassVar[int]
    KEY_FIELD_NUMBER: _ClassVar[int]
    DATA_FIELD_NUMBER: _ClassVar[int]
    index: str
    epoch: str
    create_time: str
    key: str
    data: _containers.RepeatedScalarFieldContainer[float]
    def __init__(self, index: _Optional[str] = ..., epoch: _Optional[str] = ..., create_time: _Optional[str] = ..., key: _Optional[str] = ..., data: _Optional[_Iterable[float]] = ...) -> None: ...

class LogRecord(_message.Message):
    __slots__ = ("epoch", "level", "message")
    class LogLevel(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
        __slots__ = ()
        INFO: _ClassVar[LogRecord.LogLevel]
        WARN: _ClassVar[LogRecord.LogLevel]
        ERROR: _ClassVar[LogRecord.LogLevel]
    INFO: LogRecord.LogLevel
    WARN: LogRecord.LogLevel
    ERROR: LogRecord.LogLevel
    EPOCH_FIELD_NUMBER: _ClassVar[int]
    LEVEL_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    epoch: str
    level: LogRecord.LogLevel
    message: str
    def __init__(self, epoch: _Optional[str] = ..., level: _Optional[_Union[LogRecord.LogLevel, str]] = ..., message: _Optional[str] = ...) -> None: ...
