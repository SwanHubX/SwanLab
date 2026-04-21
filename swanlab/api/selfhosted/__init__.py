"""
@author: caddiesnew
@file: project.py
@time: 2026/4/20
@description: Project 实体类 — 单个项目的查询与操作
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, cast

from swanlab.api.base import BaseEntity
from swanlab.api.typings.common import ApiResponseType
from swanlab.api.typings.selfhosted import ApiApiKeyType, ApiLicensePlanEnum, ApiSelfHostedInfoType
from swanlab.api.utils import get_properties

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class SelfHosted(BaseEntity):
    """
    表示一个 SwanLab 项目。

    支持双模式：构造时传入 data（列表迭代注入），或 data=None（按需懒加载）。
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
    ) -> None:
        super().__init__(client, web_host, api_host)

    def _ensure_data(self) -> ApiSelfHostedInfoType:
        if self._data is None:
            resp = self._get("/self_hosted/info")
            self._data = resp.data if resp.ok and resp.data else cast(ApiSelfHostedInfoType, {})
        return self._data

    @property
    def enabled(self) -> bool:
        return self._ensure_data().get("enabled", False)

    @property
    def expired(self) -> bool:
        return self._ensure_data().get("expired", False)

    @property
    def root(self) -> bool:
        return self._ensure_data().get("root", False)

    @property
    def plan(self) -> ApiLicensePlanEnum:
        return self._ensure_data().get("plan", "free")

    @property
    def seats(self) -> int:
        return self._ensure_data().get("seats", 0)

    def create_user(self, username: str, password: str) -> None:
        """
        添加用户（私有化管理员限定）
        :param username: 待创建用户名
        :param password: 待创建用户密码
        """
        data = {"users": [{"username": username, "password": password}]}
        self._post("/self_hosted/users", data=data)

    def get_users(self, page_num: int = 1, page_size: int = 20) -> ApiResponseType:
        """
        分页获取用户（管理员限定）
        :param client: 已登录的客户端实例
        :param page: 页码
        :param size: 每页大小
        """
        params = {"page": page_num, "size": page_size}
        return self._get("/self_hosted/users", params=params)
