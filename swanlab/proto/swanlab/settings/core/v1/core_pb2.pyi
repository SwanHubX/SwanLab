from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class CoreMode(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    CORE_MODE_UNSPECIFIED: _ClassVar[CoreMode]
    CORE_MODE_ONLINE: _ClassVar[CoreMode]
    CORE_MODE_LOCAL: _ClassVar[CoreMode]
    CORE_MODE_OFFLINE: _ClassVar[CoreMode]
    CORE_MODE_DISABLED: _ClassVar[CoreMode]
CORE_MODE_UNSPECIFIED: CoreMode
CORE_MODE_ONLINE: CoreMode
CORE_MODE_LOCAL: CoreMode
CORE_MODE_OFFLINE: CoreMode
CORE_MODE_DISABLED: CoreMode

class CoreSettings(_message.Message):
    __slots__ = ("run_id", "run_dir", "section_rule", "record_batch", "record_interval", "save_size", "save_split", "save_part", "save_batch", "skip_store", "api_key", "api_host", "web_host", "mode", "upload_budget_records", "upload_budget_bytes", "health_probe_interval", "health_probe_timeout", "health_degrade_threshold", "health_recover_threshold")
    RUN_ID_FIELD_NUMBER: _ClassVar[int]
    RUN_DIR_FIELD_NUMBER: _ClassVar[int]
    SECTION_RULE_FIELD_NUMBER: _ClassVar[int]
    RECORD_BATCH_FIELD_NUMBER: _ClassVar[int]
    RECORD_INTERVAL_FIELD_NUMBER: _ClassVar[int]
    SAVE_SIZE_FIELD_NUMBER: _ClassVar[int]
    SAVE_SPLIT_FIELD_NUMBER: _ClassVar[int]
    SAVE_PART_FIELD_NUMBER: _ClassVar[int]
    SAVE_BATCH_FIELD_NUMBER: _ClassVar[int]
    SKIP_STORE_FIELD_NUMBER: _ClassVar[int]
    API_KEY_FIELD_NUMBER: _ClassVar[int]
    API_HOST_FIELD_NUMBER: _ClassVar[int]
    WEB_HOST_FIELD_NUMBER: _ClassVar[int]
    MODE_FIELD_NUMBER: _ClassVar[int]
    UPLOAD_BUDGET_RECORDS_FIELD_NUMBER: _ClassVar[int]
    UPLOAD_BUDGET_BYTES_FIELD_NUMBER: _ClassVar[int]
    HEALTH_PROBE_INTERVAL_FIELD_NUMBER: _ClassVar[int]
    HEALTH_PROBE_TIMEOUT_FIELD_NUMBER: _ClassVar[int]
    HEALTH_DEGRADE_THRESHOLD_FIELD_NUMBER: _ClassVar[int]
    HEALTH_RECOVER_THRESHOLD_FIELD_NUMBER: _ClassVar[int]
    run_id: str
    run_dir: str
    section_rule: int
    record_batch: int
    record_interval: float
    save_size: int
    save_split: int
    save_part: int
    save_batch: int
    skip_store: bool
    api_key: str
    api_host: str
    web_host: str
    mode: CoreMode
    upload_budget_records: int
    upload_budget_bytes: int
    health_probe_interval: float
    health_probe_timeout: float
    health_degrade_threshold: int
    health_recover_threshold: int
    def __init__(self, run_id: _Optional[str] = ..., run_dir: _Optional[str] = ..., section_rule: _Optional[int] = ..., record_batch: _Optional[int] = ..., record_interval: _Optional[float] = ..., save_size: _Optional[int] = ..., save_split: _Optional[int] = ..., save_part: _Optional[int] = ..., save_batch: _Optional[int] = ..., skip_store: bool = ..., api_key: _Optional[str] = ..., api_host: _Optional[str] = ..., web_host: _Optional[str] = ..., mode: _Optional[_Union[CoreMode, str]] = ..., upload_budget_records: _Optional[int] = ..., upload_budget_bytes: _Optional[int] = ..., health_probe_interval: _Optional[float] = ..., health_probe_timeout: _Optional[float] = ..., health_degrade_threshold: _Optional[int] = ..., health_recover_threshold: _Optional[int] = ...) -> None: ...
