from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from typing import ClassVar as _ClassVar, Optional as _Optional

DESCRIPTOR: _descriptor.FileDescriptor

class ProbeSettings(_message.Message):
    __slots__ = ("run_id", "run_dir", "global_system_step", "hardware", "runtime", "requirements", "conda", "git", "swanlab", "monitor", "monitor_interval", "monitor_disk_dir")
    RUN_ID_FIELD_NUMBER: _ClassVar[int]
    RUN_DIR_FIELD_NUMBER: _ClassVar[int]
    GLOBAL_SYSTEM_STEP_FIELD_NUMBER: _ClassVar[int]
    HARDWARE_FIELD_NUMBER: _ClassVar[int]
    RUNTIME_FIELD_NUMBER: _ClassVar[int]
    REQUIREMENTS_FIELD_NUMBER: _ClassVar[int]
    CONDA_FIELD_NUMBER: _ClassVar[int]
    GIT_FIELD_NUMBER: _ClassVar[int]
    SWANLAB_FIELD_NUMBER: _ClassVar[int]
    MONITOR_FIELD_NUMBER: _ClassVar[int]
    MONITOR_INTERVAL_FIELD_NUMBER: _ClassVar[int]
    MONITOR_DISK_DIR_FIELD_NUMBER: _ClassVar[int]
    run_id: str
    run_dir: str
    global_system_step: int
    hardware: bool
    runtime: bool
    requirements: bool
    conda: bool
    git: bool
    swanlab: bool
    monitor: bool
    monitor_interval: int
    monitor_disk_dir: str
    def __init__(self, run_id: _Optional[str] = ..., run_dir: _Optional[str] = ..., global_system_step: _Optional[int] = ..., hardware: bool = ..., runtime: bool = ..., requirements: bool = ..., conda: bool = ..., git: bool = ..., swanlab: bool = ..., monitor: bool = ..., monitor_interval: _Optional[int] = ..., monitor_disk_dir: _Optional[str] = ...) -> None: ...
