"""
@author: Zhou QiYang
@file: user.py
@time: 2026/1/2 21:01
@description: 定义用户相关的后端API接口
"""

from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from swanlab.core_python.client import Client

from swanlab.core_python.api.type import ApiKeyType, GroupType


def create_api_key(client: "Client", *, name: str = None) -> str:
    """
    创建一个api_key，完成后返回成功信息
    :param client: 已登录的客户端实例
    :param name: api_key 的名称
    """
    if name is not None:
        data = {'name': name}
        res = client.post(f"/user/key", data=data)
    else:
        res = client.post(f"/user/key")
    return res[0]


def delete_api_key(client: "Client", *, key_id: int) -> str:
    """
    删除指定id的api_key
    :param client: 已登录的客户端实例
    :param key_id: api_key的id
    """
    res = client.delete(f"/user/key/{key_id}")
    return res[0]


def get_user_groups(client: "Client", *, username: str) -> List[GroupType]:
    """
    获取当前全部的api_key
    :param client: 已登录的客户端实例
    :param username: 用户名称
    """
    res = client.get(f"/user/{username}/groups")
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
