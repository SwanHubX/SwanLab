"""
@author: caddiesnew
@file: workspaces.py
@time: 2026/4/20
@description: Workspaces 迭代器 — 用户的工作空间列表
"""

from typing import TYPE_CHECKING, Any, Dict, Iterator, List

from .base import BaseEntity
from .workspace import Workspace

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Workspaces(BaseEntity):
    """
    用户工作空间集合的迭代器。

    用法::

        for ws in api.workspaces("username"):
            print(ws.name)
    """

    def __init__(self, client: "Client", web_host: str, api_host: str, *, username: str) -> None:
        super().__init__(client, web_host, api_host)
        self._username = username

    def _get_all_workspace_names(self) -> List[str]:
        """获取用户个人空间 + 所属团队空间名称列表。"""
        resp = self._get(f"/user/{self._username}/groups")
        group_names = [r["username"] for r in resp]
        return [self._username] + group_names

    def __iter__(self) -> Iterator[Workspace]:
        for name in self._get_all_workspace_names():
            data = self._get(f"/group/{name}")
            yield Workspace(self._client, self._web_host, self._api_host, username=name, data=data)

    def to_dict(self) -> Dict[str, Any]:
        return {"username": self._username}
