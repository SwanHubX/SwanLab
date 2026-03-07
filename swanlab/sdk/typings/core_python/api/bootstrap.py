"""
@author: cunyue
@file: auth.py
@time: 2026/3/7 18:37
@description: 鉴权相关API类型提示
"""

from typing import Optional, TypedDict


class _UserInfo(TypedDict):
    # 用户昵称
    name: Optional[str]
    # 用户名
    username: str


class LoginResponse(TypedDict):
    """登录响应"""

    # 会话ID
    sid: str
    # 过期时间，格式为 ISO 8601
    expiredAt: str
    # 用户信息
    userInfo: _UserInfo
