"""
@author: caddiesnew
@file: group.py
@time: 2026/4/22
@description: 公共查询 API 组织类型定义
"""

from typing import TypedDict

from .common import ApiGroupLiteral, ApiStatusLiteral


class ApiGroupProfileType(TypedDict):
    bio: str
    url: str
    institution: str
    school: str
    email: str
    location: str


# 在项目信息和用户信息的返回结果中，该类型的字段含义不同，注意区分
class ApiGroupType(TypedDict):
    name: str  # 组织名称 (用于user.teams)
    username: str
    comment: str
    type: ApiGroupLiteral
    status: ApiStatusLiteral
    profile: ApiGroupProfileType
