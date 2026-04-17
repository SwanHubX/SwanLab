"""
@author: cunyue
@file: __init__.py
@time: 2026/3/11 12:54
@description: SwanLab 运行时类型
"""

from typing import Literal

ResumeType = Literal["must", "allow", "never"]
"""
Run 恢复策略
"""

ModeType = Literal["disabled", "cloud", "local", "offline"]
"""
Run 运行策略
"""

FinishType = Literal["success", "crashed", "aborted"]
"""
Run 结束策略
"""


AsyncLogType = Literal["asyncio", "threading", "spawn", "fork"]
"""
异步日志上报类型
"""


WorkspaceType = Literal["TEAM", "PERSON"]
"""
工作空间类型
"""

RoleType = Literal["VISITOR", "VIEWER", "MEMBER", "OWNER"]
"""
组织成员类型
"""

IdentityType = Literal["root", "user"]
"""
Self-Hosted 用户身份类型
"""

LicensePlanType = Literal["free", "commercial"]
"""
Self-Hosted 许可证类型
"""


RunStateType = Literal["RUNNING", "FINISHED", "CRASHED", "ABORTED", "OFFLINE"]
"""
实验状态类型
"""

SidebarItemType = Literal["SCALAR", "CONFIG", "STABLE"]
"""
列过滤类型
"""

