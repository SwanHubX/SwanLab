"""
@author: Zhou Qiyang
@file: workspace
@time: 2026/1/27 20:46
@description: 工作空间相关类型
"""

from typing import TypedDict, Literal, Dict

RoleType = Literal['VISITOR', 'VIEWER', 'MEMBER', 'OWNER']


class WorkspaceType(TypedDict):
    name: str
    username: str
    profile: Dict[str, str]
    type: Literal['TEAM', 'PERSON']  # 返回的信息类型
    comment: str  # 组织或者个人的描述
    # 组织信息特有
    role: RoleType  # 组织成员的角色
