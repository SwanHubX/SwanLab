"""
@author: cunyue
@file: builder.py
@time: 2026/5/8 22:12
@description: 构建记录
"""

from typing import Optional, Union

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import (
    ChartType,
    ColumnClass,
    ColumnRecord,
    ColumnType,
    MetricColors,
    SectionType,
    YRange,
)
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord, ScalarValue
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.typings.probe_python import SystemScalar


def build_probe_column(model: SystemScalar, *, chart_index: str) -> ColumnRecord:
    y_range: Optional[YRange] = None
    if model.y_min is not None and model.y_max is not None:
        if model.y_min >= model.y_max:
            raise ValueError(f"y_min must be less than y_max, but got {model.y_min} >= {model.y_max}")
        y_range = YRange(min=model.y_min, max=model.y_max)
    metric_colors: Optional[MetricColors] = None
    if model.color is not None:
        metric_colors = MetricColors(light=model.color, dark=model.color)
    return ColumnRecord(
        column_class=ColumnClass.COLUMN_CLASS_SYSTEM,
        column_type=ColumnType.COLUMN_TYPE_SCALAR,
        column_key=helper.fmt_system_key(model.key),
        column_name=model.name,
        section_type=SectionType.SECTION_TYPE_SYSTEM,
        y_range=y_range,
        chart_index=chart_index,
        chart_name=model.chart_name,
        chart_type=ChartType.CHART_TYPE_LINE,
        metric_name=model.name,
        metric_colors=metric_colors,
    )


def build_probe_scalar(key: str, *, value: Union[int, float], step: int, timestamp: Timestamp) -> ScalarRecord:
    return ScalarRecord(
        key=helper.fmt_system_key(key),
        value=ScalarValue(number=value),
        type=ColumnType.COLUMN_TYPE_SCALAR,
        step=step,
        timestamp=timestamp,
    )
