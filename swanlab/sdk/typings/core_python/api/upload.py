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


UploadColumn = TypedDict(
    "UploadColumn",
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
        "yRange": NotRequired[Tuple[float, float]],
        "chartName": NotRequired[str],
        "chartIndex": NotRequired[str],
        "metricName": NotRequired[str],
        "metricColors": NotRequired[Tuple[str, str]],
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
# Object storage resource upload
# ============================================================


class UploadResource(TypedDict):
    """预签名 URL 文件上传描述。"""

    url: str
    source_path: str
    content_type: str


# ============================================================
# /house/metrics 信封
# ============================================================


class UploadMetricPayload(TypedDict):
    """POST /house/metrics 请求体"""

    projectId: str
    experimentId: str
    type: str
    metrics: Union[UploadScalarBatch, UploadMediaBatch, UploadLogBatch]
