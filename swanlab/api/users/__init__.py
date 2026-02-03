"""
@author: cunyue
@file: __init__.py
@time: 2026/2/3 13:30
@description: OpenApi 中的用户对象迭代器
"""

from typing import Iterator

from swanlab.api.user import User
from swanlab.core_python import Client
from swanlab.core_python.api.user import get_users


class Users:
    """
    Container for a collection of User objects.
    You can iterate over the users by for-in loop.
    """

    def __init__(self, client: Client, *, login_user: str) -> None:
        self._client = client
        self._login_user = login_user

    def __iter__(self) -> Iterator[User]:
        cur_page = 0
        while True:
            cur_page += 1
            resp = get_users(
                self._client,
                page=cur_page,
                size=20,
            )
            for u in resp['list']:
                yield User(self._client, login_user=self._login_user, username=u['username'])

            if cur_page >= resp['pages']:
                break
