"""
@author: cunyue
@file: upload.py
@time: 2026/4/23 20:01
@description: 上传相关API
"""

import sys
from typing import Any, Dict, List, Literal, Tuple, TypedDict, Union

if sys.version_info >= (3, 11):
    from typing import NotRequired, Required
else:
    from typing_extensions import NotRequired, Required

# ============================================================
# Column 列 DTO
# ============================================================


DeprecatedUploadColumn = TypedDict(
    "DeprecatedUploadColumn",
    {
        # ---- 必填 ----
        "type": Required[str],
        "key": Required[str],
        # ---- 可选 ----
        "class": NotRequired[Literal["CUSTOM", "SYSTEM"]],  # 默认为 CUSTOM
        "name": NotRequired[str],
        "error": NotRequired[Dict[str, Any]],
        "sectionName": NotRequired[str],
        "sectionType": NotRequired[Literal["PUBLIC", "SYSTEM"]],
        "yRange": NotRequired[Tuple[Union[float, None], Union[float, None]]],
        "chartName": NotRequired[str],
        "chartIndex": NotRequired[str],
        "metricName": NotRequired[str],
        "metricColors": NotRequired[Tuple[str, str]],
    },
)


DeprecatedUploadColumns = List[DeprecatedUploadColumn]
"""
列指标，多视图之前的版本
"""

UploadColumn = TypedDict(
    "UploadColumn",
    {
        "key": Required[str],
        "type": Required[str],
        "sectionName": NotRequired[str],
        "xAxis": NotRequired[str],
        "hidden": NotRequired[bool],
    },
)


class UploadColumns(TypedDict):
    """列指标 DTO，多视图之后的版本"""

    series: List[UploadColumn]


# ============================================================
# Console 日志 DTO
# ============================================================


class UploadLog(TypedDict):
    level: Literal["INFO", "ERROR"]
    epoch: int
    message: str
    create_time: str


UploadLogBatch = List[UploadLog]
"""
Console 日志指标
"""


# ============================================================
# Scalar 标量 DTO
# ============================================================


class UploadScalar(TypedDict):
    """标量指标 DTO，对应 POST /house/metrics type="scalar" 的单条 metric"""

    key: str
    index: int
    data: Union[float, Literal["nan"]]  # 字符串类型用于处理 nan 值
    create_time: NotRequired[str]


UploadScalarBatch = List[UploadScalar]
"""
标量指标
"""


# ============================================================
# Media 媒体 DTO
# ============================================================


class UploadMedia(TypedDict):
    """媒体指标 DTO，对应 POST /house/metrics type="media" 的单条 metric"""

    key: str
    index: int
    data: List[str]
    more: List[Dict[str, str]]
    create_time: NotRequired[str]


UploadMediaBatch = List[UploadMedia]
"""
媒体指标
"""


# ============================================================
# Histogram 直方图 DTO
# ============================================================


class HistogramPayload(TypedDict):
    """直方图单条 metric 的 data 字段：每桶计数与桶边界"""

    _type: NotRequired[Literal["histogram"]]  # House 侧可省略、不写进数据库的类型标记
    values: List[int]
    bins: List[float]


class UploadHistogram(TypedDict):
    """直方图指标 DTO，对应 POST /house/metrics type="histogram" 的单条 metric"""

    key: str
    epoch: Literal[0]
    index: int
    data: HistogramPayload
    create_time: NotRequired[str]


UploadHistogramBatch = List[UploadHistogram]
"""
直方图指标
"""


# ============================================================
# Object storage resource upload
# ============================================================


class UploadResource(TypedDict):
    """预签名 URL 文件上传描述。"""

    url: str
    source_path: str
    content_type: str
    size: NotRequired[int]
    tracker_key: NotRequired[str]


# ============================================================
# /house/metrics 信封
# ============================================================


class UploadMetricPayload(TypedDict):
    """POST /house/metrics 请求体"""

    projectId: str
    experimentId: str
    type: str
    metrics: Union[UploadScalarBatch, UploadMediaBatch, UploadLogBatch, UploadHistogramBatch]
