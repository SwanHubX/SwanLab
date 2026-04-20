"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: 公共查询 API 实验类型定义
"""

from typing import Dict, List, Optional, TypedDict

from swanlab.sdk.typings.run import RunStateType

from .common import ApiLabelType


class ApiExperimentUserType(TypedDict):
    username: str
    name: str


class ApiExperimentType(TypedDict):
    cuid: str
    name: str
    description: str
    labels: List[ApiLabelType]
    profile: Dict[str, object]
    show: bool
    state: RunStateType
    cluster: str
    job: str
    user: ApiExperimentUserType
    rootExpId: Optional[str]
    rootProId: Optional[str]
