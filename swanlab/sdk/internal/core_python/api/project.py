"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:00
@description: SwanLab 运行时项目API
"""

from typing import Optional

from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.core_python.api.project import ProjectType


def get_or_create_project(*, username: Optional[str], name: str, public: bool) -> ProjectType:
    """
    创建、获取项目信息，本函数还有个目标是完成多视图更新前后的兼容，breaking change 主要是在创建项目的接口上
    多视图前后的创建项目接口发生变化，原本接口不再使用，这里的流程是：
    1. 创建项目时先调用新接口，如果失败则调用旧接口
    2. 创建成功或者项目存在时再次调用获取项目信息的接口


    :param name: 项目名称
    :param username: 项目所属的用户名
    :param public: 项目是否公开
    :return: 项目信息
    """
    # 1. 参数准备
    # username 默认使用当前用户名
    username = username or client.username()
    data = {"name": name, "visibility": "PUBLIC" if public else "PRIVATE", "username": username}

    try:
        # 2. 尝试调用新接口创建项目
        # 已创建：200 ; 创建成功：201 ; 失败：4xx/5xx
        client.post(f"/projects/{username}", data=data, log_error=False)
    except ApiError:
        try:
            # 3. 如果新接口失败，尝试调用旧接口创建项目
            client.post("/project", data=data, log_error=False).data
        except ApiError as e:
            if e.response.status_code != 409:
                raise e

    # 4. 获取项目信息
    return client.get(f"/project/{username}/{name}").data
