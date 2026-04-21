"""
@author: caddiesnew
@file: project.py
@time: 2026/4/20
@description: 公共查询 API 项目类型定义
"""

from typing import Dict, List, TypedDict

from .common import ApiVisibilityEnum


class ApiProjectLabelType(TypedDict):
    name: str


class ApiProjectCountType(TypedDict):
    experiments: int
    contributors: int
    collaborators: int
    clones: int


class ApiProjectType(TypedDict):
    name: str
    username: str
    path: str
    visibility: ApiVisibilityEnum
    description: str
    group: Dict[str, str]
    projectLabels: List[ApiProjectLabelType]
    _count: ApiProjectCountType
