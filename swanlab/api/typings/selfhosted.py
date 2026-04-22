"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: 公共查询 API 私有化实例类型定义
"""

from typing import TypedDict

from .common import ApiLicensePlanLiteral


class ApiApiKeyType(TypedDict):
    id: int
    name: str
    key: str


class ApiSelfHostedInfoType(TypedDict):
    enabled: bool
    expired: bool
    root: bool
    plan: ApiLicensePlanLiteral
    seats: int
