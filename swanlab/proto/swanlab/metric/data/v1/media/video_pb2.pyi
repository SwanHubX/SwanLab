from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class VideoValue(_message.Message):
    __slots__ = ("items",)
    ITEMS_FIELD_NUMBER: _ClassVar[int]
    items: _containers.RepeatedCompositeFieldContainer[VideoItem]
    def __init__(self, items: _Optional[_Iterable[_Union[VideoItem, _Mapping]]] = ...) -> None: ...

class VideoItem(_message.Message):
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
