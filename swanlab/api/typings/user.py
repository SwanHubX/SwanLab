"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: 公共查询 API 用户类型定义
"""

from typing import TypedDict

from swanlab.sdk.typings.run import LicensePlanType


class ApiGroupType(TypedDict):
    name: str
    username: str


class ApiApiKeyType(TypedDict):
    id: int
    name: str
    key: str


class ApiSelfHostedInfoType(TypedDict):
    enabled: bool
    expired: bool
    root: bool
    plan: LicensePlanType
    seats: int
