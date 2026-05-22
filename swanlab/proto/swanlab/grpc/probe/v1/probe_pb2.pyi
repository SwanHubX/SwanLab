from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.proto.swanlab.env.v1 import metadata_pb2 as _metadata_pb2
from swanlab.proto.swanlab.settings.probe.v1 import probe_pb2 as _probe_pb2
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class DeliverProbeStartRequest(_message.Message):
    __slots__ = ("probe_settings",)
    PROBE_SETTINGS_FIELD_NUMBER: _ClassVar[int]
    probe_settings: _probe_pb2.ProbeSettings
    def __init__(self, probe_settings: _Optional[_Union[_probe_pb2.ProbeSettings, _Mapping]] = ...) -> None: ...

class GetMetadataSnapshotResponse(_message.Message):
    __slots__ = ("success", "message", "metadata")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    METADATA_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    metadata: _metadata_pb2.MetadataSnapshot
    def __init__(self, success: bool = ..., message: _Optional[str] = ..., metadata: _Optional[_Union[_metadata_pb2.MetadataSnapshot, _Mapping]] = ...) -> None: ...
