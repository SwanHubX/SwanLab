"""
@author: Zhou QiYang
@file: __init__.py.py
@time: 2026/1/10 21:44
@description: 定义用户相关的后端API接口
"""

from typing import TYPE_CHECKING, List

from swanlab.core_python.api.type import GroupType, ApiKeyType, WorkspaceType
from .self_hosted import get_self_hosted_init, create_user, get_users

if TYPE_CHECKING:
    from swanlab.core_python.client import Client


def create_api_key(client: "Client", *, name: str = None) -> None:
    """
    创建一个api_key，完成后返回成功信息
    :param client: 已登录的客户端实例
    :param name: api_key 的名称
    """
    client.post(f"/user/key", data={'name': name} if name else None)


def delete_api_key(client: "Client", *, key_id: int) -> None:
    """
    删除指定id的api_key
    :param client: 已登录的客户端实例
    :param key_id: api_key的id
    """
    client.delete(f"/user/key/{key_id}")


def get_user_groups(client: "Client", *, username: str) -> List[GroupType]:
    """
    获取用户加入的组织
    :param client: 已登录的客户端实例
    :param username: 用户名称
    """
    res = client.get(f"/user/{username}/groups")
    return res[0]


def get_workspace_info(client: "Client", *, workspace: str) -> WorkspaceType:
    """
    获取指定工作空间的信息
    :param client: 已登录的客户端实例
    :param workspace: 工作空间名称
    """
    res = client.get(f"/group/{workspace}")
    return res[0]


def get_api_keys(client: "Client") -> List[ApiKeyType]:
    """
    获取当前全部的api_key
    :param client: 已登录的客户端实例
    """
    res = client.get(f"/user/key")
    return res[0]


def get_latest_api_key(client: "Client") -> ApiKeyType:
    """
    获取最新的api_key
    :param client: 已登录的客户端实例
    """
    res = client.get(f"/user/key/latest")
    return res[0]


__all__ = [
    "create_api_key",
    "delete_api_key",
    "get_user_groups",
    "get_workspace_info",
    "get_api_keys",
    "get_latest_api_key",
    "get_self_hosted_init",
    "create_user",
    "get_users",
]
