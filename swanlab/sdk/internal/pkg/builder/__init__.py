"""
@author: cunyue
@file: __init__.py
@time: 2026/5/6 12:10
@description: 指标构建模块，输入参数，输出proto对象
本模块不负责具体的字符串长度等判断，这交给调用方处理
"""

from typing import Union

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.sdk.internal.context import RunContext

__all__ = ["build_auto_column", "build_resume_column"]


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


def build_auto_column(ctx: RunContext, data_record: Union[ScalarRecord, MediaRecord]) -> ColumnRecord:
    """
    构建一个标量列记录，此函数一般用于自动构建用户已定义的指标
    """
    # TODO: 解析 section name
    column = ColumnRecord(
        column_class=ColumnClass.COLUMN_CLASS_CUSTOM,
        column_key=data_record.key,
        column_type=data_record.type,
    )
    return column
