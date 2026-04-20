"""
@author: caddiesnew
@file: project.py
@time: 2026/4/20
@description: Project 实体类 — 单个项目的查询与操作
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional

from .base import BaseEntity
from .typings.project import ApiProjectCountType, ApiProjectType
from .utils import Label, get_properties

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Project(BaseEntity):
    """
    表示一个 SwanLab 项目。

    支持双模式：构造时传入 data（列表迭代注入），或 data=None（按需懒加载）。
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
        *,
        path: str,
        data: Optional[ApiProjectType] = None,
    ) -> None:
        super().__init__(client, web_host, api_host)
        self._path = path
        self._data = data

    def _ensure_data(self) -> ApiProjectType:
        if self._data is None:
            self._data = self._get(f"/project/{self._path}")
        return self._data

    @property
    def name(self) -> str:
        return self._ensure_data()["name"]

    @property
    def path(self) -> str:
        return self._ensure_data()["path"]

    @property
    def url(self) -> str:
        return self._build_url(f"@{self._ensure_data()['path']}")

    @property
    def description(self) -> str:
        return self._ensure_data().get("description", "")

    @property
    def visibility(self) -> str:
        return self._ensure_data().get("visibility", "PUBLIC")

    @property
    def created_at(self) -> str:
        return self._ensure_data().get("createdAt", "")

    @property
    def updated_at(self) -> str:
        return self._ensure_data().get("updatedAt", "")

    @property
    def labels(self) -> List[Label]:
        return [Label(label["name"]) for label in self._ensure_data().get("projectLabels", [])]

    @property
    def count(self) -> ApiProjectCountType:
        return self._ensure_data().get("_count", {})

    def runs(self, filters: Optional[Dict[str, object]] = None):
        """获取项目下的实验列表。"""
        from .experiments import Experiments

        return Experiments(
            self._client, self._web_host, self._api_host, path=self._ensure_data()["path"], filters=filters
        )

    def delete(self) -> None:
        """删除此项目。"""
        self._delete(f"/project/{self._ensure_data()['path']}")

    def to_dict(self) -> Dict[str, Any]:
        return get_properties(self)
