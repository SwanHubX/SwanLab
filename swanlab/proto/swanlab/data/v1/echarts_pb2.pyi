from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Optional as _Optional

DESCRIPTOR: _descriptor.FileDescriptor

class EChartsValue(_message.Message):
    __slots__ = ("path", "filename", "sha256", "size", "caption")
    PATH_FIELD_NUMBER: _ClassVar[int]
    FILENAME_FIELD_NUMBER: _ClassVar[int]
    SHA256_FIELD_NUMBER: _ClassVar[int]
    SIZE_FIELD_NUMBER: _ClassVar[int]
    CAPTION_FIELD_NUMBER: _ClassVar[int]
    path: str
    filename: str
    sha256: str
    size: int
    caption: str
    def __init__(self, path: _Optional[str] = ..., filename: _Optional[str] = ..., sha256: _Optional[str] = ..., size: _Optional[int] = ..., caption: _Optional[str] = ...) -> None: ...
