import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from swanlab.proto.swanlab.metric.column.v1 import column_pb2 as _column_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class ScalarValue(_message.Message):
    __slots__ = ("number",)
    NUMBER_FIELD_NUMBER: _ClassVar[int]
    number: float
    def __init__(self, number: _Optional[float] = ...) -> None: ...

class MediaValue(_message.Message):
    __slots__ = ("items",)
    ITEMS_FIELD_NUMBER: _ClassVar[int]
    items: _containers.RepeatedCompositeFieldContainer[MediaItem]
    def __init__(self, items: _Optional[_Iterable[_Union[MediaItem, _Mapping]]] = ...) -> None: ...

class MediaItem(_message.Message):
    __slots__ = ("filename", "sha256", "size", "caption")
    FILENAME_FIELD_NUMBER: _ClassVar[int]
    SHA256_FIELD_NUMBER: _ClassVar[int]
    SIZE_FIELD_NUMBER: _ClassVar[int]
    CAPTION_FIELD_NUMBER: _ClassVar[int]
    filename: str
    sha256: str
    size: int
    caption: str
    def __init__(self, filename: _Optional[str] = ..., sha256: _Optional[str] = ..., size: _Optional[int] = ..., caption: _Optional[str] = ...) -> None: ...

class ScalarRecord(_message.Message):
    __slots__ = ("key", "step", "type", "timestamp", "value")
    KEY_FIELD_NUMBER: _ClassVar[int]
    STEP_FIELD_NUMBER: _ClassVar[int]
    TYPE_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    VALUE_FIELD_NUMBER: _ClassVar[int]
    key: str
    step: int
    type: _column_pb2.ColumnType
    timestamp: _timestamp_pb2.Timestamp
    value: ScalarValue
    def __init__(self, key: _Optional[str] = ..., step: _Optional[int] = ..., type: _Optional[_Union[_column_pb2.ColumnType, str]] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ..., value: _Optional[_Union[ScalarValue, _Mapping]] = ...) -> None: ...

class MediaRecord(_message.Message):
    __slots__ = ("key", "step", "type", "timestamp", "value")
    KEY_FIELD_NUMBER: _ClassVar[int]
    STEP_FIELD_NUMBER: _ClassVar[int]
    TYPE_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    VALUE_FIELD_NUMBER: _ClassVar[int]
    key: str
    step: int
    type: _column_pb2.ColumnType
    timestamp: _timestamp_pb2.Timestamp
    value: MediaValue
    def __init__(self, key: _Optional[str] = ..., step: _Optional[int] = ..., type: _Optional[_Union[_column_pb2.ColumnType, str]] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ..., value: _Optional[_Union[MediaValue, _Mapping]] = ...) -> None: ...
