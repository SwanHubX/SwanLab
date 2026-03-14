import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class UpdateType(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    UPDATE_TYPE_UNSPECIFIED: _ClassVar[UpdateType]
    UPDATE_TYPE_INIT: _ClassVar[UpdateType]
    UPDATE_TYPE_PATCH: _ClassVar[UpdateType]
UPDATE_TYPE_UNSPECIFIED: UpdateType
UPDATE_TYPE_INIT: UpdateType
UPDATE_TYPE_PATCH: UpdateType

class ConfigRecord(_message.Message):
    __slots__ = ("update_type", "timestamp")
    UPDATE_TYPE_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    update_type: UpdateType
    timestamp: _timestamp_pb2.Timestamp
    def __init__(self, update_type: _Optional[_Union[UpdateType, str]] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...
