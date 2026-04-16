import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.proto.swanlab.run.v1 import run_pb2 as _run_pb2
from swanlab.proto.swanlab.metric.data.v1 import data_pb2 as _data_pb2
from swanlab.proto.swanlab.metric.column.v1 import column_pb2 as _column_pb2
from swanlab.proto.swanlab.config.v1 import config_pb2 as _config_pb2
from swanlab.proto.swanlab.system.v1 import env_pb2 as _env_pb2
from swanlab.proto.swanlab.system.v1 import console_pb2 as _console_pb2
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class Record(_message.Message):
    __slots__ = ("num", "timestamp", "start", "finish", "column", "metric", "config", "console", "metadata", "requirements", "conda")
    NUM_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    START_FIELD_NUMBER: _ClassVar[int]
    FINISH_FIELD_NUMBER: _ClassVar[int]
    COLUMN_FIELD_NUMBER: _ClassVar[int]
    METRIC_FIELD_NUMBER: _ClassVar[int]
    CONFIG_FIELD_NUMBER: _ClassVar[int]
    CONSOLE_FIELD_NUMBER: _ClassVar[int]
    METADATA_FIELD_NUMBER: _ClassVar[int]
    REQUIREMENTS_FIELD_NUMBER: _ClassVar[int]
    CONDA_FIELD_NUMBER: _ClassVar[int]
    num: int
    timestamp: _timestamp_pb2.Timestamp
    start: _run_pb2.StartRecord
    finish: _run_pb2.FinishRecord
    column: _column_pb2.ColumnRecord
    metric: _data_pb2.DataRecord
    config: _config_pb2.ConfigRecord
    console: _console_pb2.ConsoleRecord
    metadata: _env_pb2.MetadataRecord
    requirements: _env_pb2.RequirementsRecord
    conda: _env_pb2.CondaRecord
    def __init__(self, num: _Optional[int] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ..., start: _Optional[_Union[_run_pb2.StartRecord, _Mapping]] = ..., finish: _Optional[_Union[_run_pb2.FinishRecord, _Mapping]] = ..., column: _Optional[_Union[_column_pb2.ColumnRecord, _Mapping]] = ..., metric: _Optional[_Union[_data_pb2.DataRecord, _Mapping]] = ..., config: _Optional[_Union[_config_pb2.ConfigRecord, _Mapping]] = ..., console: _Optional[_Union[_console_pb2.ConsoleRecord, _Mapping]] = ..., metadata: _Optional[_Union[_env_pb2.MetadataRecord, _Mapping]] = ..., requirements: _Optional[_Union[_env_pb2.RequirementsRecord, _Mapping]] = ..., conda: _Optional[_Union[_env_pb2.CondaRecord, _Mapping]] = ...) -> None: ...
