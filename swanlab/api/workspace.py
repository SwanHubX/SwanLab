"""
@author: caddiesnew
@file: workspace.py
@time: 2026/4/20
@description: Workspace 实体类 — 工作空间的查询
"""

from typing import Any, Dict, Iterator, Optional, cast

from swanlab.api.base import ApiClientContext, BaseEntity
from swanlab.api.typings.workspace import ApiWorkspaceInfoType, ApiWorkspaceLiteral
from swanlab.api.utils import get_properties


class Workspace(BaseEntity):
    """
    表示一个 SwanLab 工作空间（个人或团队）。
    """

    def __init__(
        self,
        ctx: ApiClientContext,
        *,
        username: str,
        data: Optional[ApiWorkspaceInfoType] = None,
    ) -> None:
        super().__init__(ctx)
        self._username = username
        self._data = data

    def _ensure_data(self) -> ApiWorkspaceInfoType:
        if self._data is None:
            resp = self._get(f"/group/{self._username}")
            self._data = resp.data if resp.ok and resp.data else cast(ApiWorkspaceInfoType, {})
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
    def profile(self) -> Dict[str, str]:
        return self._ensure_data().get("profile", {})

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
    ):
        """获取工作空间下的项目列表。"""
        from swanlab.api.project import Projects

        return Projects(
            self._ctx,
            path=self.username,
            sort=sort,
            search=search,
            detail=detail,
        )

    def json(self) -> Dict[str, Any]:
        return get_properties(self)


class Workspaces(BaseEntity):
    """
    用户工作空间集合的迭代器。

    用法::

        for ws in api.workspaces("username"):
            print(ws.name)
    """

    def __init__(self, ctx: ApiClientContext, *, username: str) -> None:
        super().__init__(ctx)
        self._username = username

    def _get_all_workspace_names(self) -> list[str]:
        """获取用户个人空间 + 所属团队空间名称列表。"""
        resp = self._get(f"/user/{self._username}/groups")
        if not resp.ok:
            return [self._username]
        group_names = resp.data if isinstance(resp.data, list) else []
        return [self._username] + group_names

    def __iter__(self) -> Iterator[Workspace]:
        for name in self._get_all_workspace_names():
            resp = self._get(f"/group/{name}")
            data = resp.data if resp.ok else None
            yield Workspace(self._ctx, username=name, data=data)

    def json(self) -> Dict[str, Any]:
        return {"username": self._username}
