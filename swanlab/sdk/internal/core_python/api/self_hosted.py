"""
@author: cunyue
@file: self_hosted.py
@time: 2026/4/14 19:00
@description: SwanLab 私有化部署API
"""

from swanlab.sdk.internal.core_python import client
from swanlab.sdk.typings.core_python.api.user import SelfHostedInfoType


def get_self_hosted_init() -> SelfHostedInfoType:
    """
    获取私有化部署信息
    """
    return client.get("/self_hosted/info").data


def create_user(*, username: str, password: str) -> None:
    """
    添加用户（私有化管理员限定）
    :param username: 用户名
    :param password: 用户密码
    """
    data = {"users": [{"username": username, "password": password}]}
    client.post("/self_hosted/users", data=data)


def get_users(*, page: int = 1, size: int = 20):
    """
    分页获取用户（管理员限定）
    :param page: 页码
    :param size: 每页大小
    """
    params = {"page": page, "size": size}
    return client.get("/self_hosted/users", params=params).data
