"""
@author: cunyue
@file: user.py
@time: 2026/4/14 19:00
@description: SwanLab 运行时用户API
"""

from typing import List, Optional

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.core_python.api.user import ApiKeyType, GroupType
from swanlab.sdk.typings.core_python.api.workspace import WorkspaceInfoType


def create_api_key(*, name: Optional[str] = None) -> None:
    """
    创建一个api_key
    :param name: api_key 的名称
    """
    client.post("/user/key", data={"name": name} if name else None)


def delete_api_key(*, key_id: int) -> None:
    """
    删除指定id的api_key
    :param key_id: api_key的id
    """
    client.delete(f"/user/key/{key_id}")


def get_user_groups(*, username: str) -> List[GroupType]:
    """
    获取用户加入的组织
    :param username: 用户名称
    """
    return client.get(f"/user/{username}/groups").data


def get_workspace_info(*, path: str) -> WorkspaceInfoType:
    """
    获取指定工作空间的信息
    :param path: 工作空间名称
    """
    return client.get(f"/group/{path}").data


def get_api_keys() -> List[ApiKeyType]:
    """
    获取当前全部的api_key
    """
    return client.get("/user/key").data


def get_latest_api_key() -> ApiKeyType:
    """
    获取最新的api_key
    """
    return client.get("/user/key/latest").data
