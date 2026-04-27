"""
@author: caddiesnew
@file: project.py
@time: 2026/4/20
@description: 公共查询 API 项目类型定义
"""

from typing import Dict, List, TypedDict

from .common import ApiVisibilityLiteral


class ApiProjectLabelType(TypedDict, total=False):
    name: str
    colors: List[str]
    cuid: str


class ApiProjectCountType(TypedDict):
    experiments: int
    contributors: int
    collaborators: int
    clones: int


class ApiProjectType(TypedDict, total=False):
    cuid: str
    name: str
    username: str
    path: str
    visibility: ApiVisibilityLiteral
    description: str
    group: Dict[str, str]
    projectLabels: List[ApiProjectLabelType]
    _count: ApiProjectCountType
    createdAt: str
    updatedAt: str
    role: str
