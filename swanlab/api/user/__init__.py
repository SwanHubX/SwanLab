"""
@author: Zhou QiYang
@file: __init__.py.py
@time: 2026/1/11 16:40
@description: OpenApi 中的用户对象
"""

import re
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
from swanlab.core_python.api.user.self_hosted import create_user
from swanlab.core_python.api.user.type import SelfHostedInfoType
from swanlab.core_python.auth.providers.api_key import LoginInfo
from swanlab.core_python.client import Client

STATUS_OK = "OK"
STATUS_CREATED = "Created"


class ApiUser(ApiBase):
    def __init__(self, client: Client, login_info: LoginInfo) -> None:
        super().__init__()
        self._client = client
        self._login_info = login_info
        self._api_keys: List[ApiKeyType] = []

    @property
    def username(self) -> str:
        return self._login_info.username

    @property
    def teams(self) -> List[str]:
        resp = get_user_groups(self._client, username=self.username)
        return [r['name'] for r in resp]

    @property
    def api_keys(self) -> List[str]:
        self._api_keys = get_api_keys(self._client)
        return [r['key'] for r in self._api_keys]

    def generate_api_key(self, description: str = None) -> Optional[str]:
        api_key: Optional[ApiKeyType] = None
        res = create_api_key(self._client, name=description)
        if res == STATUS_CREATED:
            api_key = get_latest_api_key(self._client)
        return api_key['key'] if api_key else None

    def delete_api_key(self, api_key: str) -> bool:
        self._api_keys = get_api_keys(self._client)
        for key in self._api_keys:
            if key['key'] == api_key:
                res = delete_api_key(self._client, key_id=key['id'])
                if res == STATUS_OK:
                    return True
        return False


class SuperUser(ApiUser):
    def __init__(self, client: Client, login_info: LoginInfo, self_hosted: SelfHostedInfoType) -> None:
        super().__init__(client, login_info)
        self._self_hosted_info = self_hosted

    def create(self, username: str, password: str) -> bool:
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
