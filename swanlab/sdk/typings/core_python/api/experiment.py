"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API类型
"""

from typing import Optional, TypedDict


class InitExperimentType(TypedDict):
    # 实验cuid
    cuid: str
    # 实验slug(run id)
    slug: Optional[str]
    # 实验名称
    name: str
