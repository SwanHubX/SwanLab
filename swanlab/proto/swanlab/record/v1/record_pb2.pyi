import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.proto.swanlab.run.v1 import run_pb2 as _run_pb2
from swanlab.proto.swanlab.metric.data.v1 import data_pb2 as _data_pb2
from swanlab.proto.swanlab.metric.column.v1 import column_pb2 as _column_pb2
from swanlab.proto.swanlab.config.v1 import config_pb2 as _config_pb2
from swanlab.proto.swanlab.system.v1 import env_pb2 as _env_pb2
from swanlab.proto.swanlab.system.v1 import console_pb2 as _console_pb2
from swanlab.proto.swanlab.save.v1 import save_pb2 as _save_pb2
from google.protobuf.internal import containers as _containers
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Iterable as _Iterable, Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

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

class UpsertConsolesRequest(_message.Message):
    __slots__ = ("consoles",)
    CONSOLES_FIELD_NUMBER: _ClassVar[int]
    consoles: _containers.RepeatedCompositeFieldContainer[_console_pb2.ConsoleRecord]
    def __init__(self, consoles: _Optional[_Iterable[_Union[_console_pb2.ConsoleRecord, _Mapping]]] = ...) -> None: ...

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

class Record(_message.Message):
    __slots__ = ("num", "timestamp", "start", "finish", "column", "scalar", "media", "config", "console", "metadata", "requirements", "conda", "save")
    NUM_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    START_FIELD_NUMBER: _ClassVar[int]
    FINISH_FIELD_NUMBER: _ClassVar[int]
    COLUMN_FIELD_NUMBER: _ClassVar[int]
    SCALAR_FIELD_NUMBER: _ClassVar[int]
    MEDIA_FIELD_NUMBER: _ClassVar[int]
    CONFIG_FIELD_NUMBER: _ClassVar[int]
    CONSOLE_FIELD_NUMBER: _ClassVar[int]
    METADATA_FIELD_NUMBER: _ClassVar[int]
    REQUIREMENTS_FIELD_NUMBER: _ClassVar[int]
    CONDA_FIELD_NUMBER: _ClassVar[int]
    SAVE_FIELD_NUMBER: _ClassVar[int]
    num: int
    timestamp: _timestamp_pb2.Timestamp
    start: _run_pb2.StartRecord
    finish: _run_pb2.FinishRecord
    column: _column_pb2.ColumnRecord
    scalar: _data_pb2.ScalarRecord
    media: _data_pb2.MediaRecord
    config: _config_pb2.ConfigRecord
    console: _console_pb2.ConsoleRecord
    metadata: _env_pb2.MetadataRecord
    requirements: _env_pb2.RequirementsRecord
    conda: _env_pb2.CondaRecord
    save: _save_pb2.SaveRecord
    def __init__(self, num: _Optional[int] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ..., start: _Optional[_Union[_run_pb2.StartRecord, _Mapping]] = ..., finish: _Optional[_Union[_run_pb2.FinishRecord, _Mapping]] = ..., column: _Optional[_Union[_column_pb2.ColumnRecord, _Mapping]] = ..., scalar: _Optional[_Union[_data_pb2.ScalarRecord, _Mapping]] = ..., media: _Optional[_Union[_data_pb2.MediaRecord, _Mapping]] = ..., config: _Optional[_Union[_config_pb2.ConfigRecord, _Mapping]] = ..., console: _Optional[_Union[_console_pb2.ConsoleRecord, _Mapping]] = ..., metadata: _Optional[_Union[_env_pb2.MetadataRecord, _Mapping]] = ..., requirements: _Optional[_Union[_env_pb2.RequirementsRecord, _Mapping]] = ..., conda: _Optional[_Union[_env_pb2.CondaRecord, _Mapping]] = ..., save: _Optional[_Union[_save_pb2.SaveRecord, _Mapping]] = ...) -> None: ...
