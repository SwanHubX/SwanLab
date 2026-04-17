"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:00
@description: SwanLab 运行时项目API
"""

from typing import Optional, cast

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.internal.pkg.client.utils import decode_response
from swanlab.sdk.typings.core_python.api.project import InitProjectType, ProjectType, ProjResponseType


def get_project(*, username: str, name: str) -> ProjectType:
    """
    获取项目详情信息
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
        data = {"name": name, "visibility": "PUBLIC" if public else "PRIVATE", "username": username}
        return client.post("/project", data=helper.strip_none(data)).data
    except ApiError as e:
        if e.response.status_code == 409:
            # 项目已经存在，从对象中解析信息
            return cast(InitProjectType, cast(object, decode_response(e.response)))
        else:
            # 此接口为后端处理，sdk 在理论上不会出现其他错误，因此不需要处理其他错误
            raise e


def get_workspace_projects(
    *,
    path: str,
    page: int = 1,
    size: int = 20,
    sort: Optional[str] = None,
    search: Optional[str] = None,
    detail: Optional[bool] = True,
) -> ProjResponseType:
    """
    获取指定页数和条件下的项目信息
    :param path: 工作空间名称
    :param page: 页码
    :param size: 每页项目数量
    :param sort: 排序规则, 可选
    :param search: 搜索的项目名称关键字, 可选
    :param detail: 是否包含项目下实验的相关信息, 可选, 默认为true
    """
    params = {"page": page, "size": size, "sort": sort, "search": search, "detail": detail}
    return client.get(f"/project/{path}", params=helper.strip_none(params, strip_empty_str=True)).data


def delete_project(*, path: str) -> None:
    """
    删除指定项目
    :param path: 项目路径 'username/project'
    """
    client.delete(f"/project/{path}")
