"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: 公共查询 API self-hosted 类型定义
"""

from typing import TypedDict

from .common import ApiLicensePlanType


class ApiApiKeyType(TypedDict):
    id: int
    name: str
    key: str


class ApiSelfHostedInfoType(TypedDict):
    enabled: bool
    expired: bool
    root: bool
    plan: ApiLicensePlanType
    seats: int
