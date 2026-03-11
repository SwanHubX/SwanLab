import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class StreamType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    STREAM_TYPE_UNSPECIFIED: _ClassVar[StreamType]
    STREAM_TYPE_STDOUT: _ClassVar[StreamType]
    STREAM_TYPE_STDERR: _ClassVar[StreamType]
STREAM_TYPE_UNSPECIFIED: StreamType
STREAM_TYPE_STDOUT: StreamType
STREAM_TYPE_STDERR: StreamType

class ConsoleRecord(_message.Message):
    __slots__ = ("line", "stream", "timestamp")
    LINE_FIELD_NUMBER: _ClassVar[int]
    STREAM_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    line: str
    stream: StreamType
    timestamp: _timestamp_pb2.Timestamp
    def __init__(self, line: _Optional[str] = ..., stream: _Optional[_Union[StreamType, str]] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...
