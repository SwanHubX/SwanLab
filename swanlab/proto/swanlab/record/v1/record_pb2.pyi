import datetime

from google.protobuf import timestamp_pb2 as _timestamp_pb2
from google.protobuf import empty_pb2 as _empty_pb2
from swanlab.run.v1 import run_pb2 as _run_pb2
from swanlab.data.v1 import metric_pb2 as _metric_pb2
from swanlab.data.v1 import column_pb2 as _column_pb2
from swanlab.config.v1 import config_pb2 as _config_pb2
from swanlab.system.v1 import env_pb2 as _env_pb2
from swanlab.system.v1 import console_pb2 as _console_pb2
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from collections.abc import Mapping as _Mapping
from typing import ClassVar as _ClassVar, Optional as _Optional, Union as _Union

DESCRIPTOR: _descriptor.FileDescriptor

class Record(_message.Message):
    __slots__ = ("num", "timestamp", "run", "finish", "metric", "config", "column", "metadata", "requirements", "conda", "console")
    NUM_FIELD_NUMBER: _ClassVar[int]
    TIMESTAMP_FIELD_NUMBER: _ClassVar[int]
    RUN_FIELD_NUMBER: _ClassVar[int]
    FINISH_FIELD_NUMBER: _ClassVar[int]
    METRIC_FIELD_NUMBER: _ClassVar[int]
    CONFIG_FIELD_NUMBER: _ClassVar[int]
    COLUMN_FIELD_NUMBER: _ClassVar[int]
    METADATA_FIELD_NUMBER: _ClassVar[int]
    REQUIREMENTS_FIELD_NUMBER: _ClassVar[int]
    CONDA_FIELD_NUMBER: _ClassVar[int]
    CONSOLE_FIELD_NUMBER: _ClassVar[int]
    num: int
    timestamp: _timestamp_pb2.Timestamp
    run: _run_pb2.RunRecord
    finish: _run_pb2.FinishRecord
    metric: _metric_pb2.MetricRecord
    config: _config_pb2.ConfigRecord
    column: _column_pb2.ColumnRecord
    metadata: _env_pb2.MetadataRecord
    requirements: _env_pb2.RequirementsRecord
    conda: _env_pb2.CondaRecord
    console: _console_pb2.ConsoleRecord
    def __init__(self, num: _Optional[int] = ..., timestamp: _Optional[_Union[datetime.datetime, _timestamp_pb2.Timestamp, _Mapping]] = ..., run: _Optional[_Union[_run_pb2.RunRecord, _Mapping]] = ..., finish: _Optional[_Union[_run_pb2.FinishRecord, _Mapping]] = ..., metric: _Optional[_Union[_metric_pb2.MetricRecord, _Mapping]] = ..., config: _Optional[_Union[_config_pb2.ConfigRecord, _Mapping]] = ..., column: _Optional[_Union[_column_pb2.ColumnRecord, _Mapping]] = ..., metadata: _Optional[_Union[_env_pb2.MetadataRecord, _Mapping]] = ..., requirements: _Optional[_Union[_env_pb2.RequirementsRecord, _Mapping]] = ..., conda: _Optional[_Union[_env_pb2.CondaRecord, _Mapping]] = ..., console: _Optional[_Union[_console_pb2.ConsoleRecord, _Mapping]] = ...) -> None: ...
