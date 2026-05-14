from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.proto.swanlab.run.v1 import run_pb2 as _run_pb2
from swanlab.proto.swanlab.metric.data.v1 import data_pb2 as _data_pb2
from swanlab.proto.swanlab.metric.column.v1 import column_pb2 as _column_pb2
from swanlab.proto.swanlab.config.v1 import config_pb2 as _config_pb2
from swanlab.proto.swanlab.env.v1 import env_pb2 as _env_pb2
from swanlab.proto.swanlab.terminal.v1 import log_pb2 as _log_pb2
from swanlab.proto.swanlab.save.v1 import save_pb2 as _save_pb2
from swanlab.proto.swanlab.settings.core.v1 import core_pb2 as _core_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class DeliverRunStartRequest(_message.Message):
    __slots__ = ("core_settings", "start_record")
    CORE_SETTINGS_FIELD_NUMBER: _ClassVar[int]
    START_RECORD_FIELD_NUMBER: _ClassVar[int]
    core_settings: _core_pb2.CoreSettings
    start_record: _run_pb2.StartRecord
    def __init__(self, core_settings: _Optional[_Union[_core_pb2.CoreSettings, _Mapping]] = ..., start_record: _Optional[_Union[_run_pb2.StartRecord, _Mapping]] = ...) -> None: ...

class DeliverRunStartResponse(_message.Message):
    __slots__ = ("success", "message", "run", "path", "name", "global_step", "global_system_step", "new_experiment")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    RUN_FIELD_NUMBER: _ClassVar[int]
    PATH_FIELD_NUMBER: _ClassVar[int]
    NAME_FIELD_NUMBER: _ClassVar[int]
    GLOBAL_STEP_FIELD_NUMBER: _ClassVar[int]
    GLOBAL_SYSTEM_STEP_FIELD_NUMBER: _ClassVar[int]
    NEW_EXPERIMENT_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    run: _run_pb2.StartRecord
    path: str
    name: str
    global_step: int
    global_system_step: int
    new_experiment: bool
    def __init__(self, success: bool = ..., message: _Optional[str] = ..., run: _Optional[_Union[_run_pb2.StartRecord, _Mapping]] = ..., path: _Optional[str] = ..., name: _Optional[str] = ..., global_step: _Optional[int] = ..., global_system_step: _Optional[int] = ..., new_experiment: bool = ...) -> None: ...

class UpsertColumnsRequest(_message.Message):
    __slots__ = ("columns",)
    COLUMNS_FIELD_NUMBER: _ClassVar[int]
    columns: _containers.RepeatedCompositeFieldContainer[_column_pb2.ColumnRecord]
    def __init__(self, columns: _Optional[_Iterable[_Union[_column_pb2.ColumnRecord, _Mapping]]] = ...) -> None: ...

class UpsertScalarsRequest(_message.Message):
    __slots__ = ("data",)
    DATA_FIELD_NUMBER: _ClassVar[int]
    data: _containers.RepeatedCompositeFieldContainer[_data_pb2.ScalarRecord]
    def __init__(self, data: _Optional[_Iterable[_Union[_data_pb2.ScalarRecord, _Mapping]]] = ...) -> None: ...

class UpsertMediaRequest(_message.Message):
    __slots__ = ("data",)
    DATA_FIELD_NUMBER: _ClassVar[int]
    data: _containers.RepeatedCompositeFieldContainer[_data_pb2.MediaRecord]
    def __init__(self, data: _Optional[_Iterable[_Union[_data_pb2.MediaRecord, _Mapping]]] = ...) -> None: ...

class UpsertConfigsRequest(_message.Message):
    __slots__ = ("configs",)
    CONFIGS_FIELD_NUMBER: _ClassVar[int]
    configs: _containers.RepeatedCompositeFieldContainer[_config_pb2.ConfigRecord]
    def __init__(self, configs: _Optional[_Iterable[_Union[_config_pb2.ConfigRecord, _Mapping]]] = ...) -> None: ...

class UpsertLogsRequest(_message.Message):
    __slots__ = ("logs",)
    LOGS_FIELD_NUMBER: _ClassVar[int]
    logs: _containers.RepeatedCompositeFieldContainer[_log_pb2.LogRecord]
    def __init__(self, logs: _Optional[_Iterable[_Union[_log_pb2.LogRecord, _Mapping]]] = ...) -> None: ...

class UpsertMetadataRequest(_message.Message):
    __slots__ = ("metadata",)
    METADATA_FIELD_NUMBER: _ClassVar[int]
    metadata: _containers.RepeatedCompositeFieldContainer[_env_pb2.MetadataRecord]
    def __init__(self, metadata: _Optional[_Iterable[_Union[_env_pb2.MetadataRecord, _Mapping]]] = ...) -> None: ...

class UpsertRequirementsRequest(_message.Message):
    __slots__ = ("requirements",)
    REQUIREMENTS_FIELD_NUMBER: _ClassVar[int]
    requirements: _containers.RepeatedCompositeFieldContainer[_env_pb2.RequirementsRecord]
    def __init__(self, requirements: _Optional[_Iterable[_Union[_env_pb2.RequirementsRecord, _Mapping]]] = ...) -> None: ...

class UpsertCondaRequest(_message.Message):
    __slots__ = ("conda",)
    CONDA_FIELD_NUMBER: _ClassVar[int]
    conda: _containers.RepeatedCompositeFieldContainer[_env_pb2.CondaRecord]
    def __init__(self, conda: _Optional[_Iterable[_Union[_env_pb2.CondaRecord, _Mapping]]] = ...) -> None: ...

class UpsertSavesRequest(_message.Message):
    __slots__ = ("saves",)
    SAVES_FIELD_NUMBER: _ClassVar[int]
    saves: _containers.RepeatedCompositeFieldContainer[_save_pb2.SaveRecord]
    def __init__(self, saves: _Optional[_Iterable[_Union[_save_pb2.SaveRecord, _Mapping]]] = ...) -> None: ...

class DeliverRunFinishRequest(_message.Message):
    __slots__ = ("finish_record",)
    FINISH_RECORD_FIELD_NUMBER: _ClassVar[int]
    finish_record: _run_pb2.FinishRecord
    def __init__(self, finish_record: _Optional[_Union[_run_pb2.FinishRecord, _Mapping]] = ...) -> None: ...

class DeliverRunFinishResponse(_message.Message):
    __slots__ = ("success", "message")
    SUCCESS_FIELD_NUMBER: _ClassVar[int]
    MESSAGE_FIELD_NUMBER: _ClassVar[int]
    success: bool
    message: str
    def __init__(self, success: bool = ..., message: _Optional[str] = ...) -> None: ...
