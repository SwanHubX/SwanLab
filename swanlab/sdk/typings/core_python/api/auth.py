"""
@author: cunyue
@file: auth.py
@time: 2026/3/7 18:37
@description: 鉴权相关API类型提示
"""

from typing import TypedDict


class LoginResponse(TypedDict):
    """登录响应"""

    sid: str
