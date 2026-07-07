"""
@author: cunyue
@description: 项目相关将被弃用的接口
"""

import warnings

from typing_extensions import deprecated

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python import client


@deprecated("legacy projects without views will no longer be supported in v0.10.")
def get_or_create_old_project(*, data: dict):
    """
    旧接口创建、获取项目信息，主要用于兼容旧版本的项目创建接口
    :param data: 项目数据
    :return: 项目信息
    """
    warnings.warn(
        "You are accessing a legacy project that does not support views. "
        "Legacy projects will no longer be supported in v0.10, scheduled for release on 2026-10-01. "
        "Please create a new project or upgrade this project in the web app.",
        FutureWarning,
        stacklevel=2,
    )
    try:
        # 1. 尝试调用旧接口创建项目
        # 已创建：200 ; 创建成功：201 ; 失败：4xx/5xx
        client.post("/project", data=data, log_error=False)
    except ApiError as e:
        if e.response.status_code != 409:
            raise e
