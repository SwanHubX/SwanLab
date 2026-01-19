"""
@author: Zhou QiYang
@file: __init__.py
@time: 2026/1/11 16:40
@description: OpenApi 中的用户对象
"""

import re
from functools import cached_property
from typing import List, Optional

from swanlab.core_python.api.type import ApiKeyType
from swanlab.core_python.api.type.user import IdentityType
from swanlab.core_python.api.user import (
    get_user_groups,
    get_api_keys,
    create_api_key,
    get_latest_api_key,
    delete_api_key,
)
from swanlab.core_python.api.user.self_hosted import create_user
from swanlab.core_python.client import Client
from swanlab.log import swanlog


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
    def __init__(
        self, client: Client, login_user: str = None, username: str = None, identity: IdentityType = 'user'
    ) -> None:
        if login_user is None and username is None:
            raise ValueError("login_user or username are required")

        self._client = client
        self._identity = identity
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
        return [r['name'] for r in resp]

    # TODO: 管理员可以对指定用户的api_key进行操作
    @cached_property
    def api_keys(self) -> List[str]:
        """
        List of api keys the user has.
        """
        if not self.is_self:
            swanlog.warning("Getting api keys of other users has not been supported yet.")
            return []
        else:
            self._api_keys = get_api_keys(self._client)
            return [r['key'] for r in self._api_keys]

    def _refresh_api_keys(self):
        """
        Refresh the list of api keys.
        """
        del self.api_keys
        self._api_keys = get_api_keys(self._client)

    def generate_api_key(self, description: str = None) -> Optional[str]:
        """
        Generate a new api key.
        """
        if not self.is_self:
            swanlog.warning("Generating api key of other users has not been supported yet.")
            return None
        else:
            api_key: Optional[ApiKeyType] = None
            res = create_api_key(self._client, name=description)
            if res:
                api_key = get_latest_api_key(self._client)
            return api_key['key'] if api_key else None

    def delete_api_key(self, api_key: str) -> bool:
        """
        Delete an api key.
        """
        if not self.is_self:
            swanlog.warning("Deleting api key of other users has not been supported yet.")
            return False
        else:
            self._refresh_api_keys()
            for key in self._api_keys:
                if key['key'] == api_key:
                    return delete_api_key(self._client, key_id=key['id'])
            return False

    def create(self, username: str, password: str) -> bool:
        """
        Create a new user. (Only root user can create other user)
        """
        if self._identity != "root" or not self.is_self:
            swanlog.warning(f"{self._cur_username} is not allowed to create other user.")
            return False
        check_create_info(username, password)
        return create_user(self._client, username=username, password=password)


__all__ = ["User"]
