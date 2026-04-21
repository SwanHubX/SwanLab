"""
@author: caddiesnew
@file: workspace.py
@time: 2026/4/20
@description: 公共查询 API 工作空间类型定义
"""

from typing import Dict, TypedDict

from .common import ApiRoleEnum, ApiWorkspaceEnum


class ApiWorkspaceInfoType(TypedDict):
    name: str
    username: str
    profile: Dict[str, str]
    type: ApiWorkspaceEnum
    comment: str
    role: ApiRoleEnum
