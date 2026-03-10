"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:00
@description: SwanLab 运行时项目API
"""

from typing import Optional, cast

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.core_python.client.helper import decode_response
from swanlab.sdk.typings.core_python.api.project import InitProjectType, ProjectType


def get_project(*, username: str, name: str) -> ProjectType:
    """
    获取项目信息
    :param username: 项目所属的用户名
    :param name: 项目名称
    :return: 项目信息
    """
    return client.get(f"/project/{username}/{name}").data


def get_or_create_project(*, username: Optional[str], name: str, public: bool) -> InitProjectType:
    """
    创建项目，如果项目已存在，则获取项目信息
    :param name: 项目名称
    :param username: 项目所属的用户名
    :param public: 项目是否公开
    :return: 项目信息
    """
    try:
        data = {"name": name, "visibility": "PUBLIC" if public else "PRIVATE"}
        if username:
            data["username"] = username
        return client.post("/project", data=data).data
    except ApiError as e:
        if e.response.status_code == 409:
            # 项目已经存在，从对象中解析信息
            return cast(InitProjectType, cast(object, decode_response(e.response)))
        else:
            # 此接口为后端处理，sdk 在理论上不会出现其他错误，因此不需要处理其他错误
            raise e
