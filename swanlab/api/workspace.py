"""
@author: caddiesnew
@file: workspace.py
@time: 2026/4/20
@description: Workspace 实体类 — 工作空间的查询
"""

from typing import Any, Dict, Iterator, List, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.common import PaginatedQuery
from swanlab.api.typings.workspace import ApiWorkspaceLiteral, ApiWorkspaceProfileType, ApiWorkspaceType
from swanlab.api.utils import get_properties, strip_dict


class Workspace(BaseEntity):
    """
    表示一个 SwanLab 工作空间（个人或团队）。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        username: str,
        data: Optional[ApiWorkspaceType] = None,
    ) -> None:
        super().__init__(ctx)
        self._username = username
        self._data = data

    def _ensure_data(self) -> ApiWorkspaceType:
        if self._data is None:
            resp = self._get(f"/group/{self._username}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiWorkspaceType, {})
        return self._data

    @property
    def name(self) -> str:
        return self._ensure_data().get("name", "")

    @property
    def username(self) -> str:
        return self._ensure_data().get("username", "")

    @property
    def workspace_type(self) -> ApiWorkspaceLiteral:
        return self._ensure_data().get("type", "PERSON")

    @property
    def profile(self) -> Dict[str, Any]:
        return strip_dict(self._ensure_data().get("profile", {}), ApiWorkspaceProfileType)

    @property
    def comment(self) -> str:
        return self._ensure_data().get("comment", "")

    @property
    def role(self) -> str:
        return self._ensure_data().get("role", "")

    def projects(
        self,
        sort: Optional[str] = None,
        search: Optional[str] = None,
        detail: Optional[bool] = True,
        page: int = 1,
        size: int = 20,
        all: bool = False,
    ):
        from swanlab.api.project import Projects

        query = PaginatedQuery(page=page, size=size, search=search, sort=sort, all=all)
        return Projects(
            self._ctx,
            path=self.username,
            query=query,
            detail=detail,
        )

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


class Workspaces(BaseEntity):
    """
    用户工作空间集合的分页迭代器。

    用法::

        for ws in api.workspaces("username"):
            print(ws.name)
    """

    def __init__(self, ctx: ApiClientContext, *, username: str) -> None:
        super().__init__(ctx)
        self._username = username
        self._data: Optional[List[ApiWorkspaceType]] = None

    def _ensure_data(self) -> List[ApiWorkspaceType]:
        if self._data is None:
            resp = self._get(f"/user/{self._username}/groups")
            self._data = resp.data if resp.ok and resp.data else []
        assert self._data is not None
        return self._data

    def __iter__(self) -> Iterator[Workspace]:
        for item in self._ensure_data():
            yield Workspace(self._ctx, username=item["username"], data=item)

    def __len__(self) -> int:
        return len(self._ensure_data())

    def json(self) -> Dict[str, Any]:
        return {"username": self._username}
