"""
@author: caddiesnew
@file: selfhosted.py
@time: 2026/4/20
@description: SelfHosted 实体类 — 私有化部署实例的查询与管理
"""

from typing import Any, Dict, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.common import ApiResponseType
from swanlab.api.typings.selfhosted import ApiLicensePlanLiteral, ApiSelfHostedInfoType
from swanlab.api.utils import get_properties


class SelfHosted(BaseEntity):
    """
    表示一个 SwanLab 私有化部署实例。

    支持双模式：构造时传入 data，或 data=None（按需懒加载）。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        data: Optional[ApiSelfHostedInfoType] = None,
    ) -> None:
        super().__init__(ctx)
        self._data = data

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
    def plan(self) -> ApiLicensePlanLiteral:
        return self._ensure_data().get("plan", "free")

    @property
    def seats(self) -> int:
        return self._ensure_data().get("seats", 0)

    def create_user(self, username: str, password: str) -> ApiResponseType:
        """
        添加用户（私有化管理员限定）。

        :param username: 待创建用户名
        :param password: 待创建用户密码
        """
        data = {"users": [{"username": username, "password": password}]}
        return self._post("/self_hosted/users", data=data)

    def get_users(self, page: int = 1, size: int = 20) -> ApiResponseType:
        """
        分页获取用户（管理员限定）。

        :param page: 页码
        :param size: 每页大小
        """
        params = {"page": page, "size": size}
        return self._get("/self_hosted/users", params=params)

    def json(self) -> Dict[str, Any]:
        return get_properties(self)
