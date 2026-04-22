"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/21 18:40
@description: SwanLab OpenAPI 类型提示, 以 Api 前缀区分
"""

from .common import (
    ApiColumnLiteral,
    ApiIdentityLiteral,
    ApiLicensePlanLiteral,
    ApiPaginationType,
    ApiResponseType,
    ApiRoleLiteral,
    ApiRunStateLiteral,
    ApiVisibilityLiteral,
    ApiWorkspaceLiteral,
)
from .experiment import ApiExperimentLabelType, ApiExperimentType
from .project import ApiProjectCountType, ApiProjectLabelType, ApiProjectType
from .selfhosted import ApiApiKeyType, ApiSelfHostedInfoType
from .user import ApiUserProfileType, ApiUserType
from .workspace import ApiWorkspaceProfileType, ApiWorkspaceType

__all__ = [
    # Literal Definition
    "ApiColumnLiteral",
    "ApiRunStateLiteral",
    "ApiVisibilityLiteral",
    "ApiWorkspaceLiteral",
    "ApiRoleLiteral",
    "ApiIdentityLiteral",
    "ApiLicensePlanLiteral",
    # General TypedDicts
    "ApiPaginationType",
    "ApiResponseType",
    # Experiment/Run
    "ApiExperimentLabelType",
    "ApiExperimentType",
    # Project
    "ApiProjectCountType",
    "ApiProjectLabelType",
    "ApiProjectType",
    # User
    "ApiUserType",
    "ApiUserProfileType",
    # Worksapce/Group
    "ApiWorkspaceType",
    "ApiWorkspaceProfileType",
    # Misc
    "ApiApiKeyType",
    "ApiSelfHostedInfoType",
]
