from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.proto.swanlab.settings.core.v1 import core_pb2 as _core_pb2
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class DeliverSyncStartRequest(_message.Message):
    __slots__ = ("core_settings",)
    CORE_SETTINGS_FIELD_NUMBER: _ClassVar[int]
    core_settings: _core_pb2.CoreSettings
    def __init__(self, core_settings: _Optional[_Union[_core_pb2.CoreSettings, _Mapping]] = ...) -> None: ...
