"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: 公共查询 API 实验类型定义
"""

from typing import Dict, List, Optional, TypedDict

from .common import ApiLabelType, ApiRunStateType
from .user import ApiUserType


class ApiExperimentType(TypedDict):
    cuid: str
    name: str
    description: str
    labels: List[ApiLabelType]
    profile: Dict[str, object]
    show: bool
    state: ApiRunStateType
    cluster: str
    job: str
    user: ApiUserType
    rootExpId: Optional[str]
    rootProId: Optional[str]
