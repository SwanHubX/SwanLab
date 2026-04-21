"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: User 实体类 — 用户信息与 API Key 管理
"""

import re
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional

from .base import BaseEntity
from .typings.selfhosted import ApiApiKeyType, ApiSelfHostedInfoType
from .utils import get_properties

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client

_USERNAME_PATTERN = re.compile(r"^[a-zA-Z0-9_-]+$")
_PASSWORD_PATTERN = re.compile(r"^(?=.*[0-9])(?=.*[a-zA-Z]).{8,}$")


class User(BaseEntity):
    """
    表示一个 SwanLab 用户。
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
        *,
        username: str,
        login_user: str = "",
    ) -> None:
        super().__init__(client, web_host, api_host)
        self._username = username
        self._login_user = login_user or username
        self._teams: Optional[List[str]] = None
        self._api_keys_cache: Optional[List[ApiApiKeyType]] = None

    @property
    def username(self) -> str:
        return self._username

    @property
    def is_self(self) -> bool:
        return self._username == self._login_user

    @property
    def teams(self) -> List[str]:
        """用户所属的团队列表。"""
        if self._teams is None:
            resp = self._get(f"/user/{self._username}/groups")
            self._teams = resp.data if isinstance(resp.data, list) else []
        return self._teams

    @property
    def api_keys(self) -> List[str]:
        """当前用户的 API Key 列表（仅限本人）。"""
        if not self.is_self:
            return []
        if self._api_keys_cache is None:
            resp = self._get("/user/key")
            self._api_keys_cache = resp.data if resp.ok else []
        return [r["key"] for r in self._api_keys_cache or []]

    def _invalidate_api_keys_cache(self) -> None:
        self._api_keys_cache = None

    def generate_api_key(self, description: Optional[str] = None) -> Optional[str]:
        """生成新的 API Key（仅限本人）。"""
        if not self.is_self:
            return None
        self._post("/user/key", data={"name": description} if description else None)
        self._invalidate_api_keys_cache()
        resp = self._get("/user/key/latest")
        return resp.data.get("key") if resp.ok and resp.data else None

    def delete_api_key(self, api_key: str) -> bool:
        """删除指定 API Key（仅限本人）。"""
        if not self.is_self:
            return False
        self._invalidate_api_keys_cache()
        resp = self._get("/user/key")
        if not resp.ok:
            return False
        keys: List[ApiApiKeyType] = resp.data
        for key_info in keys:
            if key_info["key"] == api_key:
                self._delete(f"/user/key/{key_info['id']}")
                return True
        return False

    @staticmethod
    def _validate_create_info(username: str, password: str) -> bool:
        return bool(_USERNAME_PATTERN.match(username) and _PASSWORD_PATTERN.match(password))

    def _check_self_hosted_permission(self) -> bool:
        """检查私有化部署权限：已启用、未过期、当前用户为 root。"""
        resp = self._get("/self_hosted/info")
        if not resp.ok:
            return False
        info: ApiSelfHostedInfoType = resp.data
        return info.get("enabled", False) and not info.get("expired", True) and info.get("root", False)

    def create(self, username: str, password: str) -> bool:
        """创建新用户（仅限私有化部署 root 用户）。"""
        if (
            not self.is_self
            or not self._check_self_hosted_permission()
            or not self._validate_create_info(username, password)
        ):
            return False
        create_resp = self._post("/self_hosted/users", data={"users": [{"username": username, "password": password}]})
        return create_resp.ok

    def to_dict(self) -> Dict[str, Any]:
        return get_properties(self)


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
