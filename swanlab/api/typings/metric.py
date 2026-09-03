"""
@author: caddiesnew
@file: metric.py
@time: 2026/4/23
@description: 指标数据类型定义（用于 column 采样值）
"""

from typing import Any, List, Literal, TypedDict, Union

from swanlab.api.typings.common import ApiMetricXAxisLiteral

# ---------------------------------------------------------------------------
# Common — 通用指标类型定义
# ---------------------------------------------------------------------------


# 指标值类型（"NaN", "INF", "-INF"）
ApiMetricValueType = Union[int, float, str]


# ---------------------------------------------------------------------------
# X Axis — 查询参数与回显标记
# ---------------------------------------------------------------------------
# x_axis 查询参数：内置轴字面量，或自定义 x 列 key（任意非空字符串）
ApiMetricXAxisParam = Union[ApiMetricXAxisLiteral, str]

# 回显中 axis 的类别："step" 为内置步数轴；CUSTOM/SYSTEM 携带 key
ApiMetricXAxisKindLiteral = Literal["step", "CUSTOM", "SYSTEM"]


class ApiMetricXAxisType(TypedDict, total=False):
    """描述 ``metrics[].index`` 的实际语义，随每个 per-key metric 回显。"""

    type: ApiMetricXAxisKindLiteral
    key: str


# ---------------------------------------------------------------------------
# Column Reference — 指标列引用，标识要查询的指标列
# ---------------------------------------------------------------------------
class ApiMetricColumnRefType(TypedDict, total=False):
    projectId: str
    experimentId: str
    key: str
    rootProId: str
    rootExpId: str


# ---------------------------------------------------------------------------
# Scalar — 标量指标类型
# ---------------------------------------------------------------------------
# 使用 index 因为 x 轴可以是 step / time / relative_time / 自定义列
# 采样响应使用 data，CSV 导出记录使用 value；
# step 仅在 CSV 全量路径（all / range_query）下填充
class ApiScalarType(TypedDict, total=False):
    index: float
    data: ApiMetricValueType
    value: ApiMetricValueType
    timestamp: int
    step: int


# 组合 /metrics/scalar 和 /metrics/scalar/value 的标量序列
class ApiScalarSeriesType(ApiMetricColumnRefType, total=False):
    """标量指标序列，包含折线数据和聚合值"""

    metrics: List[ApiScalarType]
    url: str
    xAxis: ApiMetricXAxisType
    min: ApiScalarType
    max: ApiScalarType
    avg: ApiScalarType
    median: ApiScalarType
    latest: ApiScalarType


# 指标概要
# summary[run_id][key] 为下面一个 item 项
class ApiScalarSummaryItemType(TypedDict, total=False):
    step: int
    value: Any
    minMax: List[Any]
    min: Any
    max: Any
    avg: Any
    median: Any
    stdDev: Any


# ---------------------------------------------------------------------------
# Media — 媒体数据
# ---------------------------------------------------------------------------
class ApiMediaItemDataType(TypedDict, total=False):
    url: str


class ApiMediaType(TypedDict, total=False):
    index: int
    items: List[ApiMediaItemDataType]


class ApiMediaSeriesType(ApiMetricColumnRefType, total=False):
    steps: List[int]
    step: int
    metrics: List[ApiMediaType]


# ---------------------------------------------------------------------------
# Log — 日志数据
# ---------------------------------------------------------------------------
class ApiLogType(TypedDict, total=False):
    epoch: int
    level: str
    message: str
    tag: str
    timestamp: str


class ApiLogSeriesType(ApiMetricColumnRefType, total=False):
    logs: List[ApiLogType]
    count: int


# 统一数据类型定义用于类型提示
ApiMetricType = Union[ApiScalarType, ApiMediaType, ApiLogType]
