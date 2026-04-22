"""
@author: caddiesnew
@file: project.py
@time: 2026/4/20
@description: Project 实体类 — 单个项目的查询与操作
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, cast

from swanlab.api.base import BaseEntity
from swanlab.api.typings.project import ApiProjectCountType, ApiProjectLabelType, ApiProjectType
from swanlab.api.utils import get_properties

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
            resp = self._get(f"/project/{self._path}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiProjectType, {})
        return self._data

    @property
    def name(self) -> str:
        return self._ensure_data().get("name", "")

    @property
    def path(self) -> str:
        return self._ensure_data().get("path", "")

    @property
    def url(self) -> str:
        return self._build_web_url(f"@{self.path}")

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
    def labels(self) -> List[ApiProjectLabelType]:
        return [label for label in self._ensure_data().get("projectLabels", [])]

    @property
    def count(self) -> ApiProjectCountType:
        return self._ensure_data().get("_count", {})

    def runs(self, filters: Optional[Dict[str, object]] = None):
        """获取项目下的实验列表。"""
        from swanlab.api.experiment import Experiments

        return Experiments(self._client, self._web_host, self._api_host, path=self.path, filters=filters)

    def delete(self) -> bool:
        """删除此项目。"""
        resp = self._delete(f"/project/{self.path}")
        return resp.ok

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


class Projects(BaseEntity):
    """
    工作空间下项目集合的分页迭代器。

    用法::

        for project in api.projects("username"):
            print(project.name)
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
        *,
        path: str,
        sort: Optional[str] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
    ) -> None:
        super().__init__(client, web_host, api_host)
        self._path = path
        self._sort = sort
        self._search = search
        self._detail = detail

    def __iter__(self) -> Iterator[Project]:
        params = {"sort": self._sort, "search": self._search, "detail": self._detail}
        for item in self._paginate(f"/project/{self._path}", params=params):
            yield Project(
                self._client,
                self._web_host,
                self._api_host,
                path=str(item.get("path", "")),
                data=cast(ApiProjectType, item),
            )

    def json(self) -> Dict[str, Any]:
        return {"path": self._path}
