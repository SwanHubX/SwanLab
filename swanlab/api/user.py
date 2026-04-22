"""
@author: caddiesnew
@file: user.py
@time: 2026/4/20
@description: User 实体类 — 用户信息的查询
"""

from typing import Any, Dict, Optional

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.user import ApiUserProfileType
from swanlab.api.utils import get_properties, strip_dict


class User(BaseEntity):
    """
    表示一个 SwanLab 用户, 限定为通过 sdk 登录的用户。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        data: Optional[Dict[str, Any]] = None,
    ) -> None:
        super().__init__(ctx)
        self._data = data

    def _ensure_data(self) -> Dict[str, Any]:
        if self._data is None:
            resp = self._get("/user/profile")
            self._data = strip_dict(resp.data, ApiUserProfileType) if resp.ok and resp.data else {}
        return self._data

    @property
    def name(self) -> str:
        return self._ctx.name

    @property
    def username(self) -> str:
        return self._ctx.username

    @property
    def bio(self) -> str:
        return self._ensure_data().get("bio", "")

    @property
    def institution(self) -> str:
        return self._ensure_data().get("institution", "")

    @property
    def school(self) -> str:
        return self._ensure_data().get("school", "") or ""

    @property
    def email(self) -> str:
        return self._ensure_data().get("email", "") or ""

    @property
    def location(self) -> str:
        return self._ensure_data().get("location", "")

    @property
    def url(self) -> str:
        return self._ensure_data().get("url", "")

    def json(self) -> Dict[str, Any]:
        return get_properties(self)
