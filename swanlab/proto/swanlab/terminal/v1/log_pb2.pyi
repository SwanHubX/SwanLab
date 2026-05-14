import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class LogLevel(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    LOG_LEVEL_INFO: _ClassVar[LogLevel]
    LOG_LEVEL_ERROR: _ClassVar[LogLevel]
LOG_LEVEL_INFO: LogLevel
LOG_LEVEL_ERROR: LogLevel

class LogRecord(_message.Message):
    __slots__ = ("line", "level", "timestamp", "epoch")
    LINE_FIELD_NUMBER: _ClassVar[int]
    LEVEL_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    EPOCH_FIELD_NUMBER: _ClassVar[int]
    line: str
    level: LogLevel
    timestamp: _timestamp_pb2.Timestamp
    epoch: int
    def __init__(self, line: _Optional[str] = ..., level: _Optional[_Union[LogLevel, str]] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ..., epoch: _Optional[int] = ...) -> None: ...
