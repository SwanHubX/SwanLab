"""
@author: Zhou QiYang
@file: self_hosted.py
@time: 2026/1/5 17:42
@description: 私有化相关API接口
"""

from typing import TYPE_CHECKING

from swanlab.core_python.api.types import SelfHostedInfoType

if TYPE_CHECKING:
    from swanlab.core_python.client import Client


def get_self_hosted_init(client: "Client") -> SelfHostedInfoType:
    """
    获取私有化部署信息
    :param client: 已登录的客户端实例
    """
    res = client.get(f"/self_hosted/info")
    return res[0]


def create_user(client: "Client", *, username: str, password: str) -> str:
    """
    根用户添加用户
    :param client: 已登录的客户端实例
    :param username: 用户名
    :param password: 用户密码
    """
    data = {"users": [{"username": username, "password": password}]}
    res = client.post("/self_hosted/users", data=data)
    return res[0]
