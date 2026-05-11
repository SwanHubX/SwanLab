"""
@author: cunyue
@file: __init__.py
@time: 2026/5/8 21:47
@time: 2026/4/21 16:56
@description: 构建记录
"""

from typing import Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.config.v1.config_pb2 import ConfigRecord
from swanlab.proto.swanlab.env.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType, SectionType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import FinishRecord, StartRecord
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord
from swanlab.proto.swanlab.terminal.v1.log_pb2 import LogRecord
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg.counter import Counter

__all__ = [
    "build_finish_record",
    "build_start_record",
    "build_config_record",
    "build_metadata_record",
    "build_requirements_record",
    "build_conda_record",
    "build_media_record",
    "build_scalar_record",
    "build_log_record",
    "build_auto_column",
    "build_resume_column",
    "build_save_record",
]


START_RECORD_NUM = -1
"""
约定启动记录的num为-1
"""
FINISH_RECORD_NUM = -2
"""
约定结束记录的num为-2
"""
CONFIG_RECORD_NUM = -3
"""
约定配置记录的num为-3
"""
METADATA_RECORD_NUM = -4
"""
约定元信息记录的num为-4
"""
REQUIREMENTS_RECORD_NUM = -5
"""
约定依赖记录的num为-5
"""
CONDA_RECORD_NUM = -6
"""
约定conda环境记录的num为-6
"""


def _now():
    ts = Timestamp()
    ts.GetCurrentTime()
    return ts


def build_column_record(counter: Counter, column_record: ColumnRecord):
    """
    构建列记录
    """
    return Record(num=counter.inc(), column=column_record, timestamp=_now())


def build_scalar_record(counter: Counter, scalar_record: ScalarRecord):
    """
    构建数据记录
    """
    return Record(num=counter.inc(), scalar=scalar_record, timestamp=_now())


def build_media_record(counter: Counter, media_record: MediaRecord):
    """
    构建数据记录
    """
    return Record(num=counter.inc(), media=media_record, timestamp=_now())


def build_log_record(counter: Counter, epoch: Counter, log_record: LogRecord):
    """
    构建控制台记录
    """
    record = LogRecord()
    record.CopyFrom(log_record)
    record.epoch = epoch.inc()
    return Record(num=counter.inc(), log=record, timestamp=_now())


def build_start_record(start_record: StartRecord):
    """
    构建启动记录
    """
    return Record(num=START_RECORD_NUM, start=start_record, timestamp=_now())


def build_finish_record(finish_record: FinishRecord):
    """
    构建结束记录
    """
    return Record(num=FINISH_RECORD_NUM, finish=finish_record, timestamp=_now())


def build_config_record(config_record: ConfigRecord):
    """
    构建配置记录
    """
    return Record(num=CONFIG_RECORD_NUM, config=config_record, timestamp=_now())


def build_metadata_record(ts: Timestamp):
    """
    构建元信息记录
    """
    return Record(num=METADATA_RECORD_NUM, metadata=MetadataRecord(timestamp=ts), timestamp=_now())


def build_requirements_record(ts: Timestamp):
    """
    构建依赖记录
    """
    return Record(num=REQUIREMENTS_RECORD_NUM, requirements=RequirementsRecord(timestamp=ts), timestamp=_now())


def build_conda_record(ts: Timestamp):
    """
    构建conda环境记录
    """
    return Record(num=CONDA_RECORD_NUM, conda=CondaRecord(timestamp=ts), timestamp=_now())


def build_resume_column(key: str, *, media: bool = False, system: bool = False) -> ColumnRecord:
    """
    构建一个resume模式下从云端恢复的列记录
    此构建并不会恢复完整的列信息，一些不重要的，比如section name等，会被略过
    仅会恢复列的key、type、class
    :param key: 列的key
    :param media: 是否是媒体列
    :param system: 是否是系统列
    :return: ColumnRecord
    """
    if media:
        # FIXME: 暂时不知道媒体指标的类型，因此先用 COLUMN_TYPE_UNSPECIFIED 占位
        column_record = ColumnRecord(column_key=key, column_type=ColumnType.COLUMN_TYPE_UNSPECIFIED)
    else:
        column_class = ColumnClass.COLUMN_CLASS_SYSTEM if system else ColumnClass.COLUMN_CLASS_CUSTOM
        column_record = ColumnRecord(
            column_key=key, column_type=ColumnType.COLUMN_TYPE_SCALAR, column_class=column_class
        )
    return column_record


def build_auto_column(ctx: CoreContext, data_record: Union[ScalarRecord, MediaRecord]) -> ColumnRecord:
    """
    构建一个标量列记录，此函数一般用于自动构建用户已定义的指标
    """
    # Split section_name with `section_rule` setting
    section_rule_index = ctx.config.section_rule
    parts = data_record.key.split("/")
    if len(parts) >= 2:
        cut = section_rule_index % (len(parts) - 1) + 1
        section_name = "/".join(parts[:cut])
    else:
        section_name = ""

    return ColumnRecord(
        column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
        column_key=data_record.key,
        column_type=data_record.type,
        section_name=section_name,
        section_type=SectionType.SECTION_TYPE_PUBLIC,
    )


def build_save_record(counter: Counter, save_record: SaveRecord):
    """
    构建文件保存记录
    """
    return Record(num=counter.inc(), save=save_record, timestamp=_now())
