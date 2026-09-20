from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.proto.swanlab.run.v1 import run_pb2 as _run_pb2
from swanlab.proto.swanlab.metric.data.v1 import data_pb2 as _data_pb2
from swanlab.proto.swanlab.metric.column.v1 import column_pb2 as _column_pb2
from swanlab.proto.swanlab.terminal.v1 import log_pb2 as _log_pb2
from swanlab.proto.swanlab.save.v1 import save_pb2 as _save_pb2
from swanlab.proto.swanlab.operation.v1 import operation_pb2 as _operation_pb2
from swanlab.proto.swanlab.settings.core.v1 import core_pb2 as _core_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class GetCapabilitiesRequest(_message.Message):
    __slots__ = ()
    def __init__(self) -> None: ...

class GetCapabilitiesResponse(_message.Message):
    __slots__ = ("protocol_version", "core_version", "lifecycle", "supported_modes", "store", "transport", "payload", "sync", "max_receive_message_bytes")
    PROTOCOL_VERSION_FIELD_NUMBER: _ClassVar[int]
    CORE_VERSION_FIELD_NUMBER: _ClassVar[int]
    LIFECYCLE_FIELD_NUMBER: _ClassVar[int]
    SUPPORTED_MODES_FIELD_NUMBER: _ClassVar[int]
    STORE_FIELD_NUMBER: _ClassVar[int]
    TRANSPORT_FIELD_NUMBER: _ClassVar[int]
    PAYLOAD_FIELD_NUMBER: _ClassVar[int]
    SYNC_FIELD_NUMBER: _ClassVar[int]
    MAX_RECEIVE_MESSAGE_BYTES_FIELD_NUMBER: _ClassVar[int]
    protocol_version: str
    core_version: str
    lifecycle: bool
    supported_modes: _containers.RepeatedScalarFieldContainer[str]
    store: bool
    transport: bool
    payload: bool
    sync: bool
    max_receive_message_bytes: int
    def __init__(self, protocol_version: _Optional[str] = ..., core_version: _Optional[str] = ..., lifecycle: bool = ..., supported_modes: _Optional[_Iterable[str]] = ..., store: bool = ..., transport: bool = ..., payload: bool = ..., sync: bool = ..., max_receive_message_bytes: _Optional[int] = ...) -> None: ...

class TeardownServiceRequest(_message.Message):
    __slots__ = ("owner_token",)
    OWNER_TOKEN_FIELD_NUMBER: _ClassVar[int]
    owner_token: str
    def __init__(self, owner_token: _Optional[str] = ...) -> None: ...

class TeardownServiceResponse(_message.Message):
    __slots__ = ()
    def __init__(self) -> None: ...

class DeliverRunStartRequest(_message.Message):
    __slots__ = ("core_settings", "start_record")
    CORE_SETTINGS_FIELD_NUMBER: _ClassVar[int]
    START_RECORD_FIELD_NUMBER: _ClassVar[int]
    core_settings: _core_pb2.CoreSettings
    start_record: _run_pb2.StartRecord
    def __init__(self, core_settings: _Optional[_Union[_core_pb2.CoreSettings, _Mapping]] = ..., start_record: _Optional[_Union[_run_pb2.StartRecord, _Mapping]] = ...) -> None: ...

class DeliverRunStartResponse(_message.Message):
    __slots__ = ("success", "message", "run", "path", "name", "global_step", "global_system_step", "new_experiment", "run_handle")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    RUN_FIELD_NUMBER: _ClassVar[int]
    PATH_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    GLOBAL_STEP_FIELD_NUMBER: _ClassVar[int]
    GLOBAL_SYSTEM_STEP_FIELD_NUMBER: _ClassVar[int]
    NEW_EXPERIMENT_FIELD_NUMBER: _ClassVar[int]
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    run: _run_pb2.StartRecord
    path: str
    name: str
    global_step: int
    global_system_step: int
    new_experiment: bool
    run_handle: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ..., run: _Optional[_Union[_run_pb2.StartRecord, _Mapping]] = ..., path: _Optional[str] = ..., name: _Optional[str] = ..., global_step: _Optional[int] = ..., global_system_step: _Optional[int] = ..., new_experiment: bool = ..., run_handle: _Optional[str] = ...) -> None: ...

class UpsertColumnsRequest(_message.Message):
    __slots__ = ("columns", "run_handle")
    COLUMNS_FIELD_NUMBER: _ClassVar[int]
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    columns: _containers.RepeatedCompositeFieldContainer[_column_pb2.ColumnRecord]
    run_handle: str
    def __init__(self, columns: _Optional[_Iterable[_Union[_column_pb2.ColumnRecord, _Mapping]]] = ..., run_handle: _Optional[str] = ...) -> None: ...

class UpsertScalarsRequest(_message.Message):
    __slots__ = ("data", "run_handle")
    DATA_FIELD_NUMBER: _ClassVar[int]
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    data: _containers.RepeatedCompositeFieldContainer[_data_pb2.ScalarRecord]
    run_handle: str
    def __init__(self, data: _Optional[_Iterable[_Union[_data_pb2.ScalarRecord, _Mapping]]] = ..., run_handle: _Optional[str] = ...) -> None: ...

class UpsertMediaRequest(_message.Message):
    __slots__ = ("data", "run_handle")
    DATA_FIELD_NUMBER: _ClassVar[int]
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    data: _containers.RepeatedCompositeFieldContainer[_data_pb2.MediaRecord]
    run_handle: str
    def __init__(self, data: _Optional[_Iterable[_Union[_data_pb2.MediaRecord, _Mapping]]] = ..., run_handle: _Optional[str] = ...) -> None: ...

class UpsertLogsRequest(_message.Message):
    __slots__ = ("logs", "run_handle")
    LOGS_FIELD_NUMBER: _ClassVar[int]
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    logs: _containers.RepeatedCompositeFieldContainer[_log_pb2.LogRecord]
    run_handle: str
    def __init__(self, logs: _Optional[_Iterable[_Union[_log_pb2.LogRecord, _Mapping]]] = ..., run_handle: _Optional[str] = ...) -> None: ...

class UpsertSavesRequest(_message.Message):
    __slots__ = ("saves", "run_handle")
    SAVES_FIELD_NUMBER: _ClassVar[int]
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    saves: _containers.RepeatedCompositeFieldContainer[_save_pb2.SaveRecord]
    run_handle: str
    def __init__(self, saves: _Optional[_Iterable[_Union[_save_pb2.SaveRecord, _Mapping]]] = ..., run_handle: _Optional[str] = ...) -> None: ...

class DeliverRunFinishRequest(_message.Message):
    __slots__ = ("finish_record", "run_handle")
    FINISH_RECORD_FIELD_NUMBER: _ClassVar[int]
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    finish_record: _run_pb2.FinishRecord
    run_handle: str
    def __init__(self, finish_record: _Optional[_Union[_run_pb2.FinishRecord, _Mapping]] = ..., run_handle: _Optional[str] = ...) -> None: ...

class DeliverRunFinishResponse(_message.Message):
    __slots__ = ("success", "message")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ...) -> None: ...

class GetOperationStatsRequest(_message.Message):
    __slots__ = ("run_handle",)
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    run_handle: str
    def __init__(self, run_handle: _Optional[str] = ...) -> None: ...

class GetOperationStatsResponse(_message.Message):
    __slots__ = ("success", "message", "stats")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    STATS_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    stats: _operation_pb2.OperationStats
    def __init__(self, success: bool = ..., message: _Optional[str] = ..., stats: _Optional[_Union[_operation_pb2.OperationStats, _Mapping]] = ...) -> None: ...

class ConfirmRunFinishRequest(_message.Message):
    __slots__ = ("run_handle",)
    RUN_HANDLE_FIELD_NUMBER: _ClassVar[int]
    run_handle: str
    def __init__(self, run_handle: _Optional[str] = ...) -> None: ...

class ConfirmRunFinishResponse(_message.Message):
    __slots__ = ("success", "message")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ...) -> None: ...
