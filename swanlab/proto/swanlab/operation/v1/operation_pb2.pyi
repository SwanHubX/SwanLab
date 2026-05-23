from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class CoreState(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    CORE_STATE_NOT_STARTED: _ClassVar[CoreState]
    CORE_STATE_RUNNING: _ClassVar[CoreState]
    CORE_STATE_FINISHED: _ClassVar[CoreState]
CORE_STATE_NOT_STARTED: CoreState
CORE_STATE_RUNNING: CoreState
CORE_STATE_FINISHED: CoreState

class FileProgress(_message.Message):
    __slots__ = ("path", "total", "uploaded", "rate")
    PATH_FIELD_NUMBER: _ClassVar[int]
    TOTAL_FIELD_NUMBER: _ClassVar[int]
    UPLOADED_FIELD_NUMBER: _ClassVar[int]
    RATE_FIELD_NUMBER: _ClassVar[int]
    path: str
    total: int
    uploaded: int
    rate: float
    def __init__(self, path: _Optional[str] = ..., total: _Optional[int] = ..., uploaded: _Optional[int] = ..., rate: _Optional[float] = ...) -> None: ...

class OperationStats(_message.Message):
    __slots__ = ("state", "total_number", "uploaded_number", "total_records", "uploaded_records", "total_size", "uploaded_size", "rate", "files")
    STATE_FIELD_NUMBER: _ClassVar[int]
    TOTAL_NUMBER_FIELD_NUMBER: _ClassVar[int]
    UPLOADED_NUMBER_FIELD_NUMBER: _ClassVar[int]
    TOTAL_RECORDS_FIELD_NUMBER: _ClassVar[int]
    UPLOADED_RECORDS_FIELD_NUMBER: _ClassVar[int]
    TOTAL_SIZE_FIELD_NUMBER: _ClassVar[int]
    UPLOADED_SIZE_FIELD_NUMBER: _ClassVar[int]
    RATE_FIELD_NUMBER: _ClassVar[int]
    FILES_FIELD_NUMBER: _ClassVar[int]
    state: CoreState
    total_number: int
    uploaded_number: int
    total_records: int
    uploaded_records: int
    total_size: int
    uploaded_size: int
    rate: float
    files: _containers.RepeatedCompositeFieldContainer[FileProgress]
    def __init__(self, state: _Optional[_Union[CoreState, str]] = ..., total_number: _Optional[int] = ..., uploaded_number: _Optional[int] = ..., total_records: _Optional[int] = ..., uploaded_records: _Optional[int] = ..., total_size: _Optional[int] = ..., uploaded_size: _Optional[int] = ..., rate: _Optional[float] = ..., files: _Optional[_Iterable[_Union[FileProgress, _Mapping]]] = ...) -> None: ...
