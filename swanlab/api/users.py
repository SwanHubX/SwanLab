"""
@author: caddiesnew
@file: users.py
@time: 2026/4/20
@description: Users 分页迭代器 — 私有化部署管理员获取用户列表
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator

from .base import BaseEntity
from .user import User

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Users(BaseEntity):
    """
    用户集合的分页迭代器（私有化部署管理员限定）。

    用法::

        for user in api.users():
            print(user.username)
    """

    def __init__(self, client: "Client", web_host: str, api_host: str, *, login_user: str = "") -> None:
        super().__init__(client, web_host, api_host)
        self._username = login_user

    def __iter__(self) -> Iterator[User]:
        for item in self._paginate("/self_hosted/users"):
            yield User(
                self._client,
                self._web_host,
                self._api_host,
                username=item.get("username", ""),
                login_user=self._username,
            )

    def to_dict(self) -> Dict[str, Any]:
        return {}
