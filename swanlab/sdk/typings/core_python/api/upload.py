"""
@author: cunyue
@file: upload.py
@time: 2026/4/23 20:01
@description: 上传相关API
"""

import sys
from typing import Any, Dict, List, Literal, TypedDict, Union

if sys.version_info >= (3, 11):
    from typing import NotRequired, Required
else:
    from typing_extensions import NotRequired, Required

# ============================================================
# Column 列 DTO
# ============================================================


UploadColumn = TypedDict(
    "UploadColumn",
    {
        # ---- 必填 ----
        "class": Required[str],
        "type": Required[str],
        "key": Required[str],
        # ---- 可选 ----
        "name": NotRequired[str],
        "error": NotRequired[Dict[str, Any]],
        "sectionName": NotRequired[str],
        "sectionType": NotRequired[str],
        "yRange": NotRequired[Dict[str, Any]],
        "chartName": NotRequired[str],
        "chartIndex": NotRequired[str],
        "metricName": NotRequired[str],
        "metricColors": NotRequired[List[str]],
    },
)


UploadColumns = List[UploadColumn]
"""
列指标
"""


# ============================================================
# Console 日志 DTO
# ============================================================


class UploadLog(TypedDict):
    level: Literal["INFO", "ERROR"]
    epoch: int
    message: str
    create_time: str


UploadLogMetrics = List[UploadLog]
"""
Console 日志指标
"""


# ============================================================
# Scalar 标量 DTO
# ============================================================


class UploadScalar(TypedDict):
    """标量指标 DTO，对应 POST /house/metrics type="scalar" 的单条 metric"""

    key: str
    epoch: int
    index: int
    data: float
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
    epoch: int
    index: int
    data: Union[str, List[str]]
    more: NotRequired[Union[str, List[str]]]
    create_time: NotRequired[str]


UploadMediaBatch = List[UploadMedia]
"""
媒体指标
"""


# ============================================================
# /house/metrics 信封
# ============================================================


class UploadMetricPayload(TypedDict):
    """POST /house/metrics 请求体"""

    projectId: str
    experimentId: str
    type: str
    metrics: Union[UploadScalarBatch, UploadMediaBatch, UploadLogMetrics]
