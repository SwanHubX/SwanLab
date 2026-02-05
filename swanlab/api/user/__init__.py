"""
@author: Zhou QiYang
@file: __init__.py
@time: 2026/1/11 16:40
@description: OpenApi 中的用户对象
"""

import re
from functools import cached_property
from typing import List, Optional

from swanlab.api.utils import self_hosted, get_properties
from swanlab.core_python.api.type import ApiKeyType
from swanlab.core_python.api.user import (
    get_user_groups,
    get_api_keys,
    create_api_key,
    get_latest_api_key,
    delete_api_key,
    create_user,
)
from swanlab.core_python.client import Client


def check_create_info(username: str, password: str) -> bool:
    # 用户名为大小写字母、数字及-、_组成
    # 密码必须包含数字+英文且至少8位
    if not re.match(r'^[a-zA-Z0-9_-]+$', username):
        raise ValueError("Username must be alphanumeric and can contain - and _")
    if not re.match(r'^(?=.*[0-9])(?=.*[a-zA-Z]).{8,}$', password):
        raise ValueError("Password must contain at least 8 characters and include numbers and letters")
    else:
        return True


class User:
    def __init__(self, client: Client, login_user: str = None, username: str = None) -> None:
        if login_user is None and username is None:
            raise ValueError("login_user or username are required")

        self._client = client
        self._api_keys: List[ApiKeyType] = []
        self._login_user = login_user
        self._cur_username = username or self._login_user

    @property
    def username(self) -> str:
        """
        User name. (if username is not None, return username, otherwise return login_user)
        """
        return self._cur_username

    @property
    def is_self(self) -> bool:
        """
        Check if the user is the current login user.
        """
        return self._cur_username == self._login_user

    @cached_property
    def teams(self) -> List[str]:
        """
        List of teams the user belongs to.
        """
        resp = get_user_groups(self._client, username=self._cur_username)
        return [r['username'] for r in resp]

    # TODO: 管理员可以对指定用户的api_key进行操作
    @cached_property
    def api_keys(self) -> List[str]:
        """
        List of api keys the user has.
        """
        if not self.is_self:
            raise ValueError("Getting api keys of other users has not been supported yet.")
        else:
            self._api_keys = get_api_keys(self._client)
            return [r['key'] for r in self._api_keys]

    def json(self):
        """
        JSON-serializable dict of all @property values.
        """
        return get_properties(self)

    def _refresh_api_keys(self):
        """
        Refresh the list of api keys.
        """
        # Delete the cached property by removing it from instance __dict__
        self.__dict__.pop('_api_keys', None)
        self._api_keys = get_api_keys(self._client)

    def generate_api_key(self, description: str = None) -> Optional[str]:
        """
        Generate a new api key.
        """
        if not self.is_self:
            raise ValueError("Generating api key of other users has not been supported yet.")
        else:
            create_api_key(self._client, name=description)
            api_key = get_latest_api_key(self._client)
            return api_key['key'] if api_key else None

    def delete_api_key(self, api_key: str) -> bool:
        """
        Delete an api key.
        """
        if not self.is_self:
            raise ValueError("Deleting api key of other users has not been supported yet.")
        else:
            self._refresh_api_keys()
            for key in self._api_keys:
                if key['key'] == api_key:
                    delete_api_key(self._client, key_id=key['id'])
                    return True
            return False

    @self_hosted("root")
    def create(self, username: str, password: str) -> Optional[bool]:
        """
        Create a new user. (Only root user can create other user)
        """
        if not self.is_self:
            raise ValueError(f"{self._cur_username} is not allowed to create other user.")
        check_create_info(username, password)
        create_user(self._client, username=username, password=password)
        return True


__all__ = ["User"]
