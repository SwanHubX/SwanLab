"""
@author: caddiesnew
@file: workspace.py
@time: 2026/4/20
@description: Workspace 实体类 — 工作空间的查询
"""

from typing import TYPE_CHECKING, Any, Dict, Literal, Optional

from .base import BaseEntity
from .typings.workspace import ApiWorkspaceInfoType
from .utils import get_properties

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Workspace(BaseEntity):
    """
    表示一个 SwanLab 工作空间（个人或团队）。
    """

    def __init__(
        self,
        client: "Client",
        web_host: str,
        api_host: str,
        *,
        username: str,
        data: Optional[ApiWorkspaceInfoType] = None,
    ) -> None:
        super().__init__(client, web_host, api_host)
        self._username = username
        self._data = data

    def _ensure_data(self) -> ApiWorkspaceInfoType:
        if self._data is None:
            self._data = self._get(f"/group/{self._username}")
        return self._data

    @property
    def name(self) -> str:
        return self._ensure_data()["name"]

    @property
    def username(self) -> str:
        return self._ensure_data()["username"]

    @property
    def workspace_type(self) -> Literal["TEAM", "PERSON"]:
        return self._ensure_data()["type"]

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
        from .projects import Projects

        return Projects(
            self._client,
            self._web_host,
            self._api_host,
            path=self._ensure_data()["username"],
            sort=sort,
            search=search,
            detail=detail,
        )

    def to_dict(self) -> Dict[str, Any]:
        return get_properties(self)
