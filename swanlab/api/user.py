"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: User 实体类 — 用户信息的查询
"""

from typing import Any, Dict, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.user import ApiUserProfileType, ApiUserType
from swanlab.api.utils import get_properties, strip_dict


class User(BaseEntity):
    """
    表示一个 SwanLab 用户。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        data: Optional[ApiUserType] = None,
    ) -> None:
        super().__init__(ctx)
        self._data = data

    def _ensure_data(self) -> ApiUserType:
        if self._data is None:
            resp = self._get("/user/profile")
            self._data = resp.data if resp.ok and resp.data else cast(ApiUserType, {})
        return self._data

    @property
    def name(self) -> str:
        return self._ensure_data().get("name", "")

    @property
    def username(self) -> str:
        return self._ensure_data().get("username", "")

    @property
    def verified(self) -> bool:
        return self._ensure_data().get("verified", False)

    @property
    def status(self) -> str:
        return self._ensure_data().get("status", "DISABLED")

    @property
    def profile(self) -> Dict[str, Any]:
        return strip_dict(self._ensure_data().get("profile", {}), ApiUserProfileType)

    def json(self) -> Dict[str, Any]:
        return get_properties(self)
