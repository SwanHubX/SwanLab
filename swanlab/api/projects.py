"""
@author: caddiesnew
@file: projects.py
@time: 2026/4/20
@description: Projects 分页迭代器 — 工作空间下的项目列表
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator, Optional, cast

from .base import BaseEntity
from .project import Project
from .typings.project import ApiProjectType

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


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

    def to_dict(self) -> Dict[str, Any]:
        return {"path": self._path}
