"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API类型
"""

from typing import List, Optional, TypedDict


class InitExperimentType(TypedDict):
    # 实验cuid
    cuid: str
    # 实验slug(run id)
    slug: Optional[str]
    # 实验名称
    name: str


_ColumnSummary = TypedDict("_ColumnSummary", {"key": str, "step": int})


class ResumeExperimentSummaryType(TypedDict):
    """
    实验resume时的摘要信息
    """

    log: Optional[List[_ColumnSummary]]
    media: Optional[List[_ColumnSummary]]
    scalar: Optional[List[_ColumnSummary]]
