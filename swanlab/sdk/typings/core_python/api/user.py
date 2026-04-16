"""
@author: cunyue
@file: user.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时用户API类型
"""

from typing import TypedDict

from swanlab.sdk.typings.run import LicensePlanType


class GroupType(TypedDict):
    # 组织名称
    name: str
    # 组织用户名
    username: str


class ApiKeyType(TypedDict):
    # API Key ID
    id: int
    # API Key 名称
    name: str
    # API Key 值
    key: str


class SelfHostedInfoType(TypedDict):
    # 是否成功部署
    enabled: bool
    # licence是否过期
    expired: bool
    # 是否为根用户
    root: bool
    # 私有化版本（免费、商业）
    plan: LicensePlanType
    # 余剩席位
    seats: int
