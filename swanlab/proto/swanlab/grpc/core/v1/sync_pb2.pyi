from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.proto.swanlab.grpc.core.v1 import core_pb2 as _core_pb2
from swanlab.proto.swanlab.run.v1 import run_pb2 as _run_pb2
from swanlab.proto.swanlab.settings.core.v1 import core_pb2 as _core_pb2_1
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class DeliverSyncStartRequest(_message.Message):
    __slots__ = ("core_settings", "workspace", "project", "id")
    CORE_SETTINGS_FIELD_NUMBER: _ClassVar[int]
    WORKSPACE_FIELD_NUMBER: _ClassVar[int]
    PROJECT_FIELD_NUMBER: _ClassVar[int]
    ID_FIELD_NUMBER: _ClassVar[int]
    core_settings: _core_pb2_1.CoreSettings
    workspace: str
    project: str
    id: str
    def __init__(self, core_settings: _Optional[_Union[_core_pb2_1.CoreSettings, _Mapping]] = ..., workspace: _Optional[str] = ..., project: _Optional[str] = ..., id: _Optional[str] = ...) -> None: ...

class DeliverSyncStartResponse(_message.Message):
    __slots__ = ("success", "message")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ...) -> None: ...

class DeliverSyncFlushResponse(_message.Message):
    __slots__ = ("success", "message", "path")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    PATH_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    path: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ..., path: _Optional[str] = ...) -> None: ...

class ConfirmSyncFinishResponse(_message.Message):
    __slots__ = ("success", "message")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ...) -> None: ...
