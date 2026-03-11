from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class AudioValue(_message.Message):
    __slots__ = ("path", "items")
    PATH_FIELD_NUMBER: _ClassVar[int]
    ITEMS_FIELD_NUMBER: _ClassVar[int]
    path: str
    items: _containers.RepeatedCompositeFieldContainer[AudioItem]
    def __init__(self, path: _Optional[str] = ..., items: _Optional[_Iterable[_Union[AudioItem, _Mapping]]] = ...) -> None: ...

class AudioItem(_message.Message):
    __slots__ = ("filename", "sha256", "size", "format", "sample_rate", "duration_seconds", "channels", "caption")
    FILENAME_FIELD_NUMBER: _ClassVar[int]
    SHA256_FIELD_NUMBER: _ClassVar[int]
    SIZE_FIELD_NUMBER: _ClassVar[int]
    FORMAT_FIELD_NUMBER: _ClassVar[int]
    SAMPLE_RATE_FIELD_NUMBER: _ClassVar[int]
    DURATION_SECONDS_FIELD_NUMBER: _ClassVar[int]
    CHANNELS_FIELD_NUMBER: _ClassVar[int]
    CAPTION_FIELD_NUMBER: _ClassVar[int]
    filename: str
    sha256: str
    size: int
    format: str
    sample_rate: int
    duration_seconds: float
    channels: int
    caption: str
    def __init__(self, filename: _Optional[str] = ..., sha256: _Optional[str] = ..., size: _Optional[int] = ..., format: _Optional[str] = ..., sample_rate: _Optional[int] = ..., duration_seconds: _Optional[float] = ..., channels: _Optional[int] = ..., caption: _Optional[str] = ...) -> None: ...
