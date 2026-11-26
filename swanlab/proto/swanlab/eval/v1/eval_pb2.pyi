import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf.internal import enum_type_wrapper as _enum_type_wrapper
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class EvalRunState(int, metaclass=_enum_type_wrapper.EnumTypeWrapper):
    __slots__ = ()
    EVAL_RUN_STATE_RUNNING: _ClassVar[EvalRunState]
    EVAL_RUN_STATE_FINISHED: _ClassVar[EvalRunState]
    EVAL_RUN_STATE_CRASHED: _ClassVar[EvalRunState]
    EVAL_RUN_STATE_ABORTED: _ClassVar[EvalRunState]
EVAL_RUN_STATE_RUNNING: EvalRunState
EVAL_RUN_STATE_FINISHED: EvalRunState
EVAL_RUN_STATE_CRASHED: EvalRunState
EVAL_RUN_STATE_ABORTED: EvalRunState

class EvalStartRecord(_message.Message):
    __slots__ = ("project", "workspace", "cuid", "name", "tags", "link", "experiment_cuid", "step", "started_at")
    PROJECT_FIELD_NUMBER: _ClassVar[int]
    WORKSPACE_FIELD_NUMBER: _ClassVar[int]
    CUID_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    TAGS_FIELD_NUMBER: _ClassVar[int]
    LINK_FIELD_NUMBER: _ClassVar[int]
    EXPERIMENT_CUID_FIELD_NUMBER: _ClassVar[int]
    STEP_FIELD_NUMBER: _ClassVar[int]
    STARTED_AT_FIELD_NUMBER: _ClassVar[int]
    project: str
    workspace: str
    cuid: str
    name: str
    tags: _containers.RepeatedScalarFieldContainer[str]
    link: str
    experiment_cuid: str
    step: int
    started_at: _timestamp_pb2.Timestamp
    def __init__(self, project: _Optional[str] = ..., workspace: _Optional[str] = ..., cuid: _Optional[str] = ..., name: _Optional[str] = ..., tags: _Optional[_Iterable[str]] = ..., link: _Optional[str] = ..., experiment_cuid: _Optional[str] = ..., step: _Optional[int] = ..., started_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...

class EvalFinishRecord(_message.Message):
    __slots__ = ("state", "error", "finished_at")
    STATE_FIELD_NUMBER: _ClassVar[int]
    ERROR_FIELD_NUMBER: _ClassVar[int]
    FINISHED_AT_FIELD_NUMBER: _ClassVar[int]
    state: EvalRunState
    error: str
    finished_at: _timestamp_pb2.Timestamp
    def __init__(self, state: _Optional[_Union[EvalRunState, str]] = ..., error: _Optional[str] = ..., finished_at: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ...) -> None: ...
