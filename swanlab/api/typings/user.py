"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: 公共查询 API 用户类型定义
"""

from typing import TypedDict


class ApiUserType(TypedDict):
    bio: str
    institution: str
    localtion: str
    school: str
    email: str
    idc: str
    url: str
    telephone: str


class ApiGroupType(TypedDict):
    name: str
    username: str
