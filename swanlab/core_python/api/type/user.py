"""
@author: Zhou QiYang
@file: user.py
@time: 2026/1/10 21:46
@description: 用户相关后端API接口类型
"""

from typing import Literal, TypedDict


# ------------------------------------- 通用类型 -------------------------------------
IdentityType = Literal['user', 'root']


# 在项目信息和用户信息的返回结果中，该类型的字段含义不同，注意区分
class GroupType(TypedDict):
    name: str  # 组织名称 (用于user.teams)
    username: str


# ------------------------------------- 后端返回信息 -------------------------------------
class ApiKeyType(TypedDict):
    id: int
    name: str
    createdAt: str
    key: str


# 私有化部署信息
class SelfHostedInfoType(TypedDict):
    enabled: bool  # 是否成功部署
    expired: bool  # licence是否过期
    root: bool  # 是否为根用户
    plan: Literal["free", "commercial"]  # 私有化版本（免费、商业）
    seats: int  # 余剩席位
