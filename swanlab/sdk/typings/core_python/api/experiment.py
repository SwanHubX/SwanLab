"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API类型
"""

from typing import List, Optional, TypedDict

from typing_extensions import NotRequired

# createdAt 参数的合法区间（10 位 Unix Timestamp）
CREATED_AT_MIN = 1_000_000_000
CREATED_AT_MAX = 9_999_999_999


class InitExperimentType(TypedDict):
    # 实验cuid
    cuid: str
    # 实验slug(run id)
    slug: Optional[str]
    # 实验名称
    name: str
    # 实验创建时间，ISO 8601 字符串，resume 时作为 House 查询的 createdAt 下界
    createdAt: NotRequired[str]


_ColumnSummary = TypedDict("_ColumnSummary", {"key": str, "step": int})


class ResumeExperimentSummaryType(TypedDict):
    """
    实验resume时的摘要信息
    """

    log: Optional[List[_ColumnSummary]]
    media: Optional[List[_ColumnSummary]]
    scalar: Optional[List[_ColumnSummary]]
