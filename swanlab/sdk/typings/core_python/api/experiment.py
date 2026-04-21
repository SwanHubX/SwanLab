"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API类型
"""

from typing import List, Tuple, TypedDict

from swanlab.sdk.typings.run import RunStateType

from .common import LabelType


class InitExperimentType(TypedDict):
    # 实验cuid
    cuid: str


class ExperimentType(TypedDict):
    # 实验名称
    name: str
    # 实验CUID, 唯一标识符
    cuid: str
    # 实验描述
    description: str
    # 实验标签列表
    labels: List[LabelType]
    # 实验状态
    state: RunStateType
    # 实验组
    cluster: str
    # 任务类型
    job: str
    # 创建时间
    created_at: str
    # 实验颜色
    colors: Tuple[str, str]
