"""
@author: caddiesnew
@file: group.py
@time: 2026/4/20
@description: Group 实体类 — 组织信息的查询
"""

from typing import Any, Dict, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.group import ApiGroupProfileType, ApiGroupType
from swanlab.api.utils import get_properties, strip_dict


class Group(BaseEntity):
    """
    表示一个 SwanLab 组织。
    """

    def __init__(self, ctx: ApiClientContext, username: str, data: Optional[ApiGroupType] = None) -> None:
        super().__init__(ctx)
        self._username = username
        self._data = data

    def _ensure_data(self) -> ApiGroupType:
        if self._data is None:
            resp = self._get(f"/group/{self._username}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiGroupType, {})
        return self._data

    @property
    def name(self) -> str:
        return self._ensure_data().get("name", "")

    @property
    def username(self) -> str:
        return self._ensure_data().get("username", "")

    @property
    def comment(self) -> str:
        return self._ensure_data().get("comment", "")

    @property
    def group_type(self) -> str:
        return self._ensure_data().get("type", "TEAM")

    @property
    def status(self) -> str:
        return self._ensure_data().get("status", "ACTIVE")

    @property
    def profile(self) -> Dict[str, Any]:
        return strip_dict(self._ensure_data().get("profile", {}), ApiGroupProfileType)

    def json(self) -> Dict[str, Any]:
        return get_properties(self)
