"""
@author: cunyue
@file: __init__.py
@time: 2026/5/6 12:10
@description: 指标构建模块，输入参数，输出proto对象
本模块不负责具体的字符串长度等判断，这交给调用方处理
"""

from typing import Union

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaRecord, ScalarRecord
from swanlab.sdk.internal.context import RunContext

__all__ = ["build_auto_column"]


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
