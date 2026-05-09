from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class SavePolicy(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    SAVE_POLICY_NOW: _ClassVar[SavePolicy]
    SAVE_POLICY_END: _ClassVar[SavePolicy]
    SAVE_POLICY_LIVE: _ClassVar[SavePolicy]
SAVE_POLICY_NOW: SavePolicy
SAVE_POLICY_END: SavePolicy
SAVE_POLICY_LIVE: SavePolicy

class SaveRecord(_message.Message):
    __slots__ = ("name", "source_path", "target_path", "policy")
    NAME_FIELD_NUMBER: _ClassVar[int]
    SOURCE_PATH_FIELD_NUMBER: _ClassVar[int]
    TARGET_PATH_FIELD_NUMBER: _ClassVar[int]
    POLICY_FIELD_NUMBER: _ClassVar[int]
    name: str
    source_path: str
    target_path: str
    policy: SavePolicy
    def __init__(self, name: _Optional[str] = ..., source_path: _Optional[str] = ..., target_path: _Optional[str] = ..., policy: _Optional[_Union[SavePolicy, str]] = ...) -> None: ...
