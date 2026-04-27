"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/21 18:40
@description: SwanLab OpenAPI 类型提示, 以 Api 前缀区分
"""

from .column import ApiColumnCsvExportType, ApiColumnErrorType, ApiColumnType
from .common import (
    ApiIdentityLiteral,
    ApiLicensePlanLiteral,
    ApiMetricAllTypeLiteral,
    ApiMetricColumnTypeLiteral,
    ApiMetricLogLevelLiteral,
    ApiMetricXAxisLiteral,
    ApiPaginationType,
    ApiResponseType,
    ApiRoleLiteral,
    ApiRunStateLiteral,
    ApiSidebarLiteral,
    ApiVisibilityLiteral,
    ApiWorkspaceLiteral,
)
from .experiment import ApiExperimentLabelType, ApiExperimentType
from .metric import (
    ApiLogSeriesType,
    ApiMediaSeriesType,
    ApiScalarSeriesType,
)
from .project import ApiProjectCountType, ApiProjectLabelType, ApiProjectType
from .selfhosted import ApiApiKeyType, ApiSelfHostedInfoType
from .user import ApiUserProfileType, ApiUserType
from .workspace import ApiWorkspaceProfileType, ApiWorkspaceType

__all__ = [
    # Literal Definition
    "ApiSidebarLiteral",
    "ApiRunStateLiteral",
    "ApiVisibilityLiteral",
    "ApiWorkspaceLiteral",
    "ApiRoleLiteral",
    "ApiIdentityLiteral",
    "ApiLicensePlanLiteral",
    "ApiMetricLogLevelLiteral",
    "ApiMetricAllTypeLiteral",
    "ApiMetricColumnTypeLiteral",
    "ApiMetricXAxisLiteral",
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
    # Column
    "ApiColumnErrorType",
    "ApiColumnType",
    "ApiColumnCsvExportType",
    # Metric
    "ApiLogSeriesType",
    "ApiMediaSeriesType",
    "ApiScalarSeriesType",
]
