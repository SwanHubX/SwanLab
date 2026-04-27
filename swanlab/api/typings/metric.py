"""
@author: caddiesnew
@file: metric.py
@time: 2026/4/23
@description: 指标数据类型定义（用于 column 采样值）
"""

from typing import Any, List, TypedDict, Union

# ---------------------------------------------------------------------------
# Common — 通用指标类型定义
# ---------------------------------------------------------------------------


# 指标值类型（"NaN", "INF", "-INF"）
ApiMetricValueType = Union[int, float, str]


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
class ApiScalarType(TypedDict, total=False):
    index: float
    data: ApiMetricValueType
    timestamp: int


# 组合 /metrics/scalar 和 /metrics/scalar/value 的标量序列
class ApiScalarSeriesType(ApiMetricColumnRefType, total=False):
    """标量指标序列，包含折线数据和聚合值"""

    metrics: List[ApiScalarType]
    url: str
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
