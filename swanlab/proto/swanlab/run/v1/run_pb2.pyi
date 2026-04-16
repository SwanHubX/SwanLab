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
    RUN_STATE_RUNNING: _ClassVar[RunState]
    RUN_STATE_FINISHED: _ClassVar[RunState]
    RUN_STATE_CRASHED: _ClassVar[RunState]
    RUN_STATE_ABORTED: _ClassVar[RunState]

class ResumeMode(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    RESUME_MODE_NEVER: _ClassVar[ResumeMode]
    RESUME_MODE_ALLOW: _ClassVar[ResumeMode]
    RESUME_MODE_MUST: _ClassVar[ResumeMode]
RUN_STATE_RUNNING: RunState
RUN_STATE_FINISHED: RunState
RUN_STATE_CRASHED: RunState
RUN_STATE_ABORTED: RunState
RESUME_MODE_NEVER: ResumeMode
RESUME_MODE_ALLOW: ResumeMode
RESUME_MODE_MUST: ResumeMode

class StartRequest(_message.Message):
    __slots__ = ("project", "workspace", "name", "description", "tags", "group", "job_type", "id", "resume", "started_at")
    PROJECT_FIELD_NUMBER: _ClassVar[int]
    WORKSPACE_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    TAGS_FIELD_NUMBER: _ClassVar[int]
    GROUP_FIELD_NUMBER: _ClassVar[int]
    JOB_TYPE_FIELD_NUMBER: _ClassVar[int]
    ID_FIELD_NUMBER: _ClassVar[int]
    RESUME_FIELD_NUMBER: _ClassVar[int]
    STARTED_AT_FIELD_NUMBER: _ClassVar[int]
    project: str
    workspace: str
    name: str
    description: str
    tags: _containers.RepeatedScalarFieldContainer[str]
    group: str
    job_type: str
    id: str
    resume: ResumeMode
    started_at: _timestamp_pb2.Timestamp
    def __init__(self, project: _Optional[str] = ..., workspace: _Optional[str] = ..., name: _Optional[str] = ..., description: _Optional[str] = ..., tags: _Optional[_Iterable[str]] = ..., group: _Optional[str] = ..., job_type: _Optional[str] = ..., id: _Optional[str] = ..., resume: _Optional[_Union[ResumeMode, str]] = ..., started_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...

class StartResponse(_message.Message):
    __slots__ = ("success", "message", "color")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    color: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ..., color: _Optional[str] = ...) -> None: ...

class StartRecord(_message.Message):
    __slots__ = ("project", "workspace", "name", "description", "tags", "group", "job_type", "color", "id", "resume", "started_at")
    PROJECT_FIELD_NUMBER: _ClassVar[int]
    WORKSPACE_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    DESCRIPTION_FIELD_NUMBER: _ClassVar[int]
    TAGS_FIELD_NUMBER: _ClassVar[int]
    GROUP_FIELD_NUMBER: _ClassVar[int]
    JOB_TYPE_FIELD_NUMBER: _ClassVar[int]
    COLOR_FIELD_NUMBER: _ClassVar[int]
    ID_FIELD_NUMBER: _ClassVar[int]
    RESUME_FIELD_NUMBER: _ClassVar[int]
    STARTED_AT_FIELD_NUMBER: _ClassVar[int]
    project: str
    workspace: str
    name: str
    description: str
    tags: _containers.RepeatedScalarFieldContainer[str]
    group: str
    job_type: str
    color: str
    id: str
    resume: ResumeMode
    started_at: _timestamp_pb2.Timestamp
    def __init__(self, project: _Optional[str] = ..., workspace: _Optional[str] = ..., name: _Optional[str] = ..., description: _Optional[str] = ..., tags: _Optional[_Iterable[str]] = ..., group: _Optional[str] = ..., job_type: _Optional[str] = ..., color: _Optional[str] = ..., id: _Optional[str] = ..., resume: _Optional[_Union[ResumeMode, str]] = ..., started_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...

class FinishRequest(_message.Message):
    __slots__ = ("state", "error", "finished_at")
    STATE_FIELD_NUMBER: _ClassVar[int]
    ERROR_FIELD_NUMBER: _ClassVar[int]
    FINISHED_AT_FIELD_NUMBER: _ClassVar[int]
    state: RunState
    error: str
    finished_at: _timestamp_pb2.Timestamp
    def __init__(self, state: _Optional[_Union[RunState, str]] = ..., error: _Optional[str] = ..., finished_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...

class FinishResponse(_message.Message):
    __slots__ = ("success", "message")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ...) -> None: ...

class FinishRecord(_message.Message):
    __slots__ = ("state", "error", "finished_at")
    STATE_FIELD_NUMBER: _ClassVar[int]
    ERROR_FIELD_NUMBER: _ClassVar[int]
    FINISHED_AT_FIELD_NUMBER: _ClassVar[int]
    state: RunState
    error: str
    finished_at: _timestamp_pb2.Timestamp
    def __init__(self, state: _Optional[_Union[RunState, str]] = ..., error: _Optional[str] = ..., finished_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...
