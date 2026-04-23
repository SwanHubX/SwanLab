"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: 公共查询 API 实验类型定义
"""

from typing import Any, Dict, List, Optional, TypedDict

from .common import ApiExperimentTypeLiteral, ApiRunStateLiteral
from .user import ApiUserType


class ApiExperimentLabelType(TypedDict):
    name: str


# 实验配置
class ApiExperimentProfileType(TypedDict):
    config: Dict[str, Any]
    metadata: Dict[str, Any]
    requirements: str
    conda: str


class ApiExperimentType(TypedDict):
    cuid: str
    name: str
    type: ApiExperimentTypeLiteral
    description: str
    labels: List[ApiExperimentLabelType]
    profile: ApiExperimentProfileType
    show: bool
    state: ApiRunStateLiteral
    cluster: str
    job: str
    user: ApiUserType
    rootExpId: Optional[str]
    rootProId: Optional[str]
