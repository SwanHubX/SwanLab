"""
@author: cunyue
@file: project.py
@time: 2026/3/10 18:00
@description: SwanLab 运行时项目API
"""

from typing import Optional

from swanlab.deprecated.project import get_or_create_old_project
from swanlab.exceptions import ApiError
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.core_python.api.project import ProjectType


def get_or_create_project(*, username: Optional[str], name: str, public: bool) -> ProjectType:
    """
    创建、获取项目信息，本函数还有个目标是完成多视图更新前后的兼容，breaking change 主要是在创建项目的接口上
    多视图前后的创建项目接口发生变化，原本接口不再使用，这里的流程是：
    1. 先调用获取项目信息的接口，项目已存在则直接返回（避免只写权限用户触发创建接口的 403）
    2. 项目不存在（404）时调用新接口创建项目，如果新接口不可用则调用旧接口
    3. 创建成功后再次调用获取项目信息的接口


    :param name: 项目名称
    :param username: 项目所属的用户名
    :param public: 项目是否公开
    :return: 项目信息
    """
    # 1. 参数准备
    # username 默认使用当前用户名
    username = username or client.username()
    data = {"name": name, "visibility": "PUBLIC" if public else "PRIVATE"}

    # 2. 先请求接口获取项目信息，如果404再尝试创建项目
    # 因为此时项目可能是已经存在的，但是当前用户只有写入权限，此时调用创建项目的接口返回403，属于非预期行为
    try:
        return client.get(f"/project/{username}/{name}", log_error=False).data
    except ApiError as e:
        if e.response.status_code != 404:
            raise e
    try:
        # 3. 尝试调用新接口创建项目
        # 已创建：200 ; 创建成功：201 ; 失败：4xx/5xx
        client.post(f"/projects/{username}", data=data)
    except ApiError as e:
        # 4. 如果新接口不存在，则调用旧接口创建项目
        if e.response.status_code == 404:
            get_or_create_old_project(data={**data, "username": username})
        else:
            raise e

    # 5. 获取项目信息
    return client.get(f"/project/{username}/{name}").data
