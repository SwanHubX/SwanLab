"""
@author: caddiesnew
@file: workspace.py
@time: 2026/4/20
@description: 公共查询 API 工作空间类型定义
"""

from typing import TypedDict

from .common import ApiRoleLiteral, ApiWorkspaceLiteral


# 工作空间即 Group 组织
class ApiWorkspaceProfileType(TypedDict):
    bio: str
    url: str
    institution: str
    school: str
    email: str
    location: str


class ApiWorkspaceType(TypedDict):
    username: str
    name: str
    type: ApiWorkspaceLiteral
    comment: str
    role: ApiRoleLiteral
    profile: ApiWorkspaceProfileType
