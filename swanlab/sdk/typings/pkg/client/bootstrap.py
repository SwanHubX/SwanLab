"""
@author: cunyue
@file: bootstrap.py
@time: 2026/3/7 18:37
@description: 鉴权相关API类型提示
"""

from typing import TypedDict


class UserProfile(TypedDict):
    """`GET /api/auth/verify` 返回的用户档案"""

    # 用户ID
    uid: int
    # 头像
    avatar: str
    # 用户名
    username: str
    # 用户昵称
    name: str
    # 创建时间，格式为 ISO 8601
    createdAt: str
    # 是否已验证
    verified: bool
