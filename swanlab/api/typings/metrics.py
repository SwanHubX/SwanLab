"""
@author: caddiesnew
@file: metrics.py
@time: 2026/4/23
@description: 指标数据类型定义（用于 column 采样值获取）
"""

from typing import Any, Dict, List, TypedDict


# ---------------------------------------------------------------------------
# Scalar — 标量 item 数据
# ---------------------------------------------------------------------------
# 一个 key 下单个 scalar 的数据值包装
class ApiScalarType(TypedDict, total=False):
    step: int
    data: float
    timestamp: int


# 一个 key 下批量 scalar 的响应值包装
# 需要请求两次
class ApiScalarListType(TypedDict, total=False):
    min: ApiScalarType
    max: ApiScalarType
    avg: ApiScalarType
    median: ApiScalarType
    latest: ApiScalarType
    metrics: List[ApiScalarType]


# 指标概要
# summary[run_id][key] 为下面一个 item 项
class ApiSummaryItemType(TypedDict, total=False):
    step: int
    value: Any
    minMax: List[Any]
    min: Any
    max: Any
    avg: Any
    median: Any
    stdDev: Any


# ---------------------------------------------------------------------------
# Media — 媒体 item 数据
# ---------------------------------------------------------------------------
class ApiMediaType(TypedDict, total=False):
    # 项目路径: proj_id/run_id 拼接而成
    prefix: str
    data: List[str]
    more: List[Dict[str, Any]]


# ---------------------------------------------------------------------------
# Log — 日志 item 数据
# ---------------------------------------------------------------------------
class ApiLogType(TypedDict, total=False):
    epoch: int
    level: str
    message: str
    tag: str
    timestamp: str
