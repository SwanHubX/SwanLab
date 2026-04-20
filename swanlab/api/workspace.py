"""
@author: Nexisato
@file: workspace.py
@time: 2026/4/20
@description: Workspace 实体类 — 工作空间的查询
"""

from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional

from swanlab.sdk.typings.core_python.api.workspace import WorkspaceInfoType

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Workspace:
    """
    表示一个 SwanLab 工作空间（个人或团队）。

    通过独立 Client 实例与云端交互，支持查询工作空间详情及项目列表。
    """

    def __init__(self, client: "Client", *, username: str, web_host: str) -> None:
        self._client: "Client" = client
        self._username: str = username
        self._web_host: str = web_host
        self._data: Optional[WorkspaceInfoType] = None

    def _fetch(self) -> WorkspaceInfoType:
        """从云端拉取工作空间详情，缓存到 _data。"""
        raise NotImplementedError("Workspace._fetch")

    # ------------------------------------------------------------------
    #  属性占位
    # ------------------------------------------------------------------

    @property
    def name(self) -> str:
        """工作空间显示名称。"""
        raise NotImplementedError

    @property
    def username(self) -> str:
        """工作空间标识（用户名或团队名）。"""
        raise NotImplementedError

    @property
    def workspace_type(self) -> Literal["TEAM", "PERSON"]:
        """工作空间类型。"""
        raise NotImplementedError

    # ------------------------------------------------------------------
    #  操作占位
    # ------------------------------------------------------------------

    def projects(self) -> List[Any]:
        """获取工作空间下的项目列表。"""
        raise NotImplementedError("Workspace.projects")

    def json(self) -> Dict[str, Any]:
        """返回所有属性的 JSON 可序列化字典。"""
        raise NotImplementedError("Workspace.json")


__all__ = ["Workspace"]
