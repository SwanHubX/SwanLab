from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Optional as _Optional

DESCRIPTOR: _descriptor.FileDescriptor

class CoreSettings(_message.Message):
    __slots__ = ("run_id", "run_dir", "section_rule", "record_batch", "record_interval", "save_size", "save_split", "save_part", "save_batch")
    RUN_ID_FIELD_NUMBER: _ClassVar[int]
    RUN_DIR_FIELD_NUMBER: _ClassVar[int]
    SECTION_RULE_FIELD_NUMBER: _ClassVar[int]
    RECORD_BATCH_FIELD_NUMBER: _ClassVar[int]
    RECORD_INTERVAL_FIELD_NUMBER: _ClassVar[int]
    SAVE_SIZE_FIELD_NUMBER: _ClassVar[int]
    SAVE_SPLIT_FIELD_NUMBER: _ClassVar[int]
    SAVE_PART_FIELD_NUMBER: _ClassVar[int]
    SAVE_BATCH_FIELD_NUMBER: _ClassVar[int]
    run_id: str
    run_dir: str
    section_rule: int
    record_batch: int
    record_interval: float
    save_size: int
    save_split: int
    save_part: int
    save_batch: int
    def __init__(self, run_id: _Optional[str] = ..., run_dir: _Optional[str] = ..., section_rule: _Optional[int] = ..., record_batch: _Optional[int] = ..., record_interval: _Optional[float] = ..., save_size: _Optional[int] = ..., save_split: _Optional[int] = ..., save_part: _Optional[int] = ..., save_batch: _Optional[int] = ...) -> None: ...
