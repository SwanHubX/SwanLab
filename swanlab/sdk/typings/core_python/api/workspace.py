"""
@author: cunyue
@file: workspace.py
@time: 2026/3/10 19:02
@description: SwanLab 运行时工作空间API类型
"""

from typing import Dict, TypedDict

from swanlab.sdk.typings.run import RoleType, WorkspaceType


class WorkspaceInfoType(TypedDict):
    # 工作空间名称
    name: str
    # 工作空间用户名
    username: str
    # 工作空间配置信息
    profile: Dict[str, str]
    # 工作空间类型
    type: WorkspaceType
    # 组织或者个人的描述
    comment: str
    # 组织成员的角色（组织信息特有）
    role: RoleType
