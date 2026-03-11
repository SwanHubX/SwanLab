import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class RunState(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    RUN_STATE_UNSPECIFIED: _ClassVar[RunState]
    RUN_STATE_RUNNING: _ClassVar[RunState]
    RUN_STATE_FINISHED: _ClassVar[RunState]
    RUN_STATE_CRASHED: _ClassVar[RunState]
    RUN_STATE_STOPPED: _ClassVar[RunState]

class ResumeMode(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    RESUME_MODE_UNSPECIFIED: _ClassVar[ResumeMode]
    RESUME_MODE_NEVER: _ClassVar[ResumeMode]
    RESUME_MODE_ALLOW: _ClassVar[ResumeMode]
    RESUME_MODE_MUST: _ClassVar[ResumeMode]
RUN_STATE_UNSPECIFIED: RunState
RUN_STATE_RUNNING: RunState
RUN_STATE_FINISHED: RunState
RUN_STATE_CRASHED: RunState
RUN_STATE_STOPPED: RunState
RESUME_MODE_UNSPECIFIED: ResumeMode
RESUME_MODE_NEVER: ResumeMode
RESUME_MODE_ALLOW: ResumeMode
RESUME_MODE_MUST: ResumeMode

class RunRecord(_message.Message):
    __slots__ = ("run_id", "experiment_name", "project_name", "workspace", "description", "tags", "group", "job_type", "color", "resume", "started_at")
    RUN_ID_FIELD_NUMBER: _ClassVar[int]
    EXPERIMENT_NAME_FIELD_NUMBER: _ClassVar[int]
    PROJECT_NAME_FIELD_NUMBER: _ClassVar[int]
    WORKSPACE_FIELD_NUMBER: _ClassVar[int]
    DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    TAGS_FIELD_NUMBER: _ClassVar[int]
    GROUP_FIELD_NUMBER: _ClassVar[int]
    JOB_TYPE_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    RESUME_FIELD_NUMBER: _ClassVar[int]
    STARTED_AT_FIELD_NUMBER: _ClassVar[int]
    run_id: str
    experiment_name: str
    project_name: str
    workspace: str
    description: str
    tags: _containers.RepeatedScalarFieldContainer[str]
    group: str
    job_type: str
    color: str
    resume: ResumeMode
    started_at: _timestamp_pb2.Timestamp
    def __init__(self, run_id: _Optional[str] = ..., experiment_name: _Optional[str] = ..., project_name: _Optional[str] = ..., workspace: _Optional[str] = ..., description: _Optional[str] = ..., tags: _Optional[_Iterable[str]] = ..., group: _Optional[str] = ..., job_type: _Optional[str] = ..., color: _Optional[str] = ..., resume: _Optional[_Union[ResumeMode, str]] = ..., started_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...

class FinishRecord(_message.Message):
    __slots__ = ("state", "exit_code", "error", "finished_at")
    STATE_FIELD_NUMBER: _ClassVar[int]
    EXIT_CODE_FIELD_NUMBER: _ClassVar[int]
    ERROR_FIELD_NUMBER: _ClassVar[int]
    FINISHED_AT_FIELD_NUMBER: _ClassVar[int]
    state: RunState
    exit_code: int
    error: str
    finished_at: _timestamp_pb2.Timestamp
    def __init__(self, state: _Optional[_Union[RunState, str]] = ..., exit_code: _Optional[int] = ..., error: _Optional[str] = ..., finished_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...
