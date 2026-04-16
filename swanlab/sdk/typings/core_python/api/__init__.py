"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 18:40
@description: SwanLab API类型提示

所有后端响应类型命名以 Type 结尾
"""

from .experiment import RunType
from .project import InitProjectType, ProjectLabelType, ProjectType, ProjResponseType
from .user import ApiKeyType, GroupType, SelfHostedInfoType
from .workspace import WorkspaceInfoType

__all__ = [
    "RunType",
    "ProjectType",
    "InitProjectType",
    "ProjResponseType",
    "ProjectLabelType",
    "GroupType",
    "ApiKeyType",
    "SelfHostedInfoType",
    "WorkspaceInfoType",
]
