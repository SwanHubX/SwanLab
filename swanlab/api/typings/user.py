"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: 公共查询 API 用户类型定义
"""

from typing import TypedDict


class ApiUserProfileType(TypedDict):
    # 简介
    bio: str
    # 个人链接
    url: str
    # 机构
    institution: str
    # 学校
    school: str
    # 邮箱
    email: str
    # 地址
    location: str


class ApiUserType(TypedDict):
    name: str
    username: str
