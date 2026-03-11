from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Optional as _Optional

DESCRIPTOR: _descriptor.FileDescriptor

class ScalarValue(_message.Message):
    __slots__ = ("number", "string_val", "bool_val", "is_nan", "is_inf")
    NUMBER_FIELD_NUMBER: _ClassVar[int]
    STRING_VAL_FIELD_NUMBER: _ClassVar[int]
    BOOL_VAL_FIELD_NUMBER: _ClassVar[int]
    IS_NAN_FIELD_NUMBER: _ClassVar[int]
    IS_INF_FIELD_NUMBER: _ClassVar[int]
    number: float
    string_val: str
    bool_val: bool
    is_nan: bool
    is_inf: bool
    def __init__(self, number: _Optional[float] = ..., string_val: _Optional[str] = ..., bool_val: bool = ..., is_nan: bool = ..., is_inf: bool = ...) -> None: ...
