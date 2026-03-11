from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class TextValue(_message.Message):
    __slots__ = ("path", "items")
    PATH_FIELD_NUMBER: _ClassVar[int]
    ITEMS_FIELD_NUMBER: _ClassVar[int]
    path: str
    items: _containers.RepeatedCompositeFieldContainer[TextItem]
    def __init__(self, path: _Optional[str] = ..., items: _Optional[_Iterable[_Union[TextItem, _Mapping]]] = ...) -> None: ...

class TextItem(_message.Message):
    __slots__ = ("filename", "sha256", "caption")
    FILENAME_FIELD_NUMBER: _ClassVar[int]
    SHA256_FIELD_NUMBER: _ClassVar[int]
    CAPTION_FIELD_NUMBER: _ClassVar[int]
    filename: str
    sha256: str
    caption: str
    def __init__(self, filename: _Optional[str] = ..., sha256: _Optional[str] = ..., caption: _Optional[str] = ...) -> None: ...
