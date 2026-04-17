"""
@author: cunyue
@file: experiment.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时实验API类型
"""

from typing import Dict, List, Optional, TypedDict

from swanlab.sdk.typings.run import RunStateType


class _ExperimentLabelType(TypedDict):
    # 标签名称
    name: str


class _UserType(TypedDict):
    # 用户名
    username: str
    # 用户显示名称
    name: str


class RunType(TypedDict):
    # 实验CUID, 唯一标识符
    cuid: str
    # 实验名称
    name: str
    # 实验描述
    description: str
    # 实验标签列表
    labels: List[_ExperimentLabelType]
    # 实验配置和摘要信息，包含 'config' 和 'scalar'
    profile: Dict[str, Dict[str, object]]
    # 是否显示
    show: bool
    # 实验状态
    state: RunStateType
    # 实验组
    cluster: str
    # 任务类型
    job: str
    # 实验所属用户
    user: _UserType
    # 祖宗实验对应的实验cuid，如果为克隆实验则必传
    rootExpId: Optional[str]
    # 祖宗实验对应的项目cuid，如果为克隆实验则必传
    rootProId: Optional[str]
