"""
@author: Zhou QiYang
@file: __init__.py.py
@time: 2026/1/11 16:40
@description: OpenApi 中的用户对象
"""

import re
from functools import cached_property
from typing import List, Optional

from swanlab.api.base import ApiBase
from swanlab.core_python.api.user import (
    get_user_groups,
    get_api_keys,
    create_api_key,
    get_latest_api_key,
    delete_api_key,
    ApiKeyType,
)
from swanlab.core_python.api.user.self_hosted import create_user, get_self_hosted_init
from swanlab.core_python.client import Client
from swanlab.error import ApiError
from swanlab.log import swanlog

STATUS_OK = "OK"
STATUS_CREATED = "Created"


class User(ApiBase):
    def __init__(self, client: Client, login_user: str = None, username: str = None) -> None:
        if login_user is None and username is None:
            raise ValueError("login_user or username are required")

        super().__init__()
        self._client = client
        self._api_keys: List[ApiKeyType] = []
        self._cur_username = username or login_user
        self._is_other_user = username is not None and username != login_user

        # 尝试获取私有化服务信息，如果不是私有化服务，则会报错退出，因为指定user功能仅供私有化用户使用
        try:
            self._self_hosted_info = get_self_hosted_init(self._client)
        except ApiError:
            swanlog.warning("You haven't launched a swanlab self-hosted instance. Some usages are not available.")
            self._self_hosted_info = None

        if self._self_hosted_info is not None and self._self_hosted_info["plan"] == 'commercial':
            self._identity = 'root' if self._self_hosted_info['root'] else 'user'

        if self._self_hosted_info is not None:
            if not self._self_hosted_info["enabled"]:
                swanlog.warning("SwanLab self-hosted instance hasn't been ready yet.")
            if self._self_hosted_info["expired"]:
                swanlog.warning("SwanLab self-hosted instance has expired.")

    @property
    def username(self) -> str:
        return self._cur_username

    @cached_property
    def teams(self) -> List[str]:
        resp = get_user_groups(self._client, username=self._cur_username)
        return [r['name'] for r in resp]

    # TODO: 管理员可以对指定用户的api_key进行操作
    @cached_property
    def api_keys(self) -> List[str]:
        if self._is_other_user:
            swanlog.warning("Getting api keys of other users has not been supported yet.")
            return []
        else:
            self._api_keys = get_api_keys(self._client)
            return [r['key'] for r in self._api_keys]

    def _refresh_api_keys(self):
        del self.api_keys
        self._api_keys = get_api_keys(self._client)

    def generate_api_key(self, description: str = None) -> Optional[str]:
        if self._is_other_user:
            swanlog.warning("Generating api key of other users has not been supported yet.")
            return None
        else:
            api_key: Optional[ApiKeyType] = None
            res = create_api_key(self._client, name=description)
            if res == STATUS_CREATED:
                api_key = get_latest_api_key(self._client)
            return api_key['key'] if api_key else None

    def delete_api_key(self, api_key: str) -> bool:
        if self._is_other_user:
            swanlog.warning("Deleting api key of other users has not been supported yet.")
            return False
        else:
            self._refresh_api_keys()
            for key in self._api_keys:
                if key['key'] == api_key:
                    res = delete_api_key(self._client, key_id=key['id'])
                    if res == STATUS_OK:
                        return True
            return False

    def create(self, username: str, password: str) -> bool:
        if self._identity != "root" or self._is_other_user:
            swanlog.warning(f"{self._cur_username} is not allowed to create other user.")
            return False
        # 用户名为大小写字母、数字及-、_组成
        # 密码必须包含数字+英文且至少8位
        if not re.match(r'^[a-zA-Z0-9_-]+$', username):
            raise ValueError("Username must be alphanumeric and can contain - and _")
        if not re.match(r'^(?=.*[0-9])(?=.*[a-zA-Z]).{8,}$', password):
            raise ValueError("Password must contain at least 8 characters and include numbers and letters")
        resp = create_user(self._client, username=username, password=password)
        if resp == STATUS_CREATED:
            return True
        else:
            raise False


__all__ = ["User"]
