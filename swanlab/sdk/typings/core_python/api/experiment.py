"""
@author: cunyue
@file: experiment.py
@time: 2026/4/18 17:43
@description: SwanLab 运行时实验API类型
"""

from typing import TypedDict


class InitExperimentType(TypedDict):
    # 实验cuid
    cuid: str
