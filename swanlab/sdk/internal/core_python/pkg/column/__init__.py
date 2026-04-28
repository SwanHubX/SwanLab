"""
@author: cunyue
@file: __init__.py
@time: 2026/4/25 16:17
@description: 列处理模块
"""

from typing import Literal, Optional, cast

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.sdk.internal.pkg import adapter, safe
from swanlab.sdk.typings.core_python.api.upload import UploadColumn


@safe.decorator(message="Failed to encode column record")
def encode(record: ColumnRecord) -> Optional[UploadColumn]:
    """
    将列记录编码为后端所需的格式（DTO）
    """
    column: UploadColumn = {"key": record.column_key, "type": adapter.column[record.column_type]}
    # class: 是否为系统列
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["class"] = cast(Literal["SYSTEM"], "SYSTEM")
    # name: 列的展示名称
    if record.column_name:
        column["name"] = record.column_name
    # section name: section 名称
    if record.section_name:
        column["sectionName"] = record.section_name
    # section type: section 类型，目前仅处理系统列
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["sectionType"] = cast(Literal["SYSTEM"], "SYSTEM")
    # yRange: 数值列的 y 轴范围
    if record.column_type == ColumnType.COLUMN_TYPE_SCALAR:
        if record.HasField("y_range"):
            column["yRange"] = (record.y_range.min, record.y_range.max)
    # chartName: 图表名称
    if record.chart_name:
        column["chartName"] = record.chart_name
    # chartIndex: 图表索引
    if record.chart_index:
        column["chartIndex"] = record.chart_index
    # metricName: 指标名称
    if record.metric_name:
        column["metricName"] = record.metric_name
    # metricColors: 指标颜色
    if record.HasField("metric_colors"):
        column["metricColors"] = (record.metric_colors.light, record.metric_colors.dark)
    return column


@safe.decorator(message="Failed to decode column record")
def decode(data: bytes):
    """
    将后端返回的数据解码为列记录
    """
    raise NotImplementedError("Waiting for resume feature")
