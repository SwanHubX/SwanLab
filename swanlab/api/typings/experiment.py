"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: 公共查询 API 实验类型定义
"""

from typing import Dict, List, Optional, TypedDict

from .common import ApiRunStateLiteral
from .user import ApiUserType


class ApiExperimentLabelType(TypedDict):
    name: str


class ApiExperimentType(TypedDict):
    cuid: str
    name: str
    description: str
    labels: List[ApiExperimentLabelType]
    profile: Dict[str, object]
    show: bool
    state: ApiRunStateLiteral
    cluster: str
    job: str
    user: ApiUserType
    rootExpId: Optional[str]
    rootProId: Optional[str]
