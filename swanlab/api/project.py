"""
@author: Nexisato
@file: project.py
@time: 2026/4/20
@description: Project 实体类 — 单个项目的查询与操作
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional

from swanlab.sdk.typings.core_python.api.project import ProjectType

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Project:
    """
    表示一个 SwanLab 项目。

    通过独立 Client 实例与云端交互，支持查询项目详情、获取实验列表、删除等操作。
    """

    def __init__(self, client: "Client", *, path: str, web_host: str) -> None:
        self._client: "Client" = client
        self._path: str = path  # 'username/project-name'
        self._web_host: str = web_host
        self._data: Optional[ProjectType] = None

    def _fetch(self) -> ProjectType:
        """从云端拉取项目详情，缓存到 _data。"""
        raise NotImplementedError("Project._fetch")

    # ------------------------------------------------------------------
    #  属性占位
    # ------------------------------------------------------------------

    @property
    def name(self) -> str:
        """项目名称。"""
        raise NotImplementedError

    @property
    def path(self) -> str:
        """项目路径（username/project-name）。"""
        raise NotImplementedError

    @property
    def url(self) -> str:
        """项目在 Web 面板中的访问 URL。"""
        raise NotImplementedError

    # ------------------------------------------------------------------
    #  操作占位
    # ------------------------------------------------------------------

    def runs(self, filters: Optional[Dict[str, Any]] = None) -> List[Any]:
        """获取项目下的实验列表。"""
        raise NotImplementedError("Project.runs")

    def delete(self) -> None:
        """删除此项目。"""
        raise NotImplementedError("Project.delete")

    def json(self) -> Dict[str, Any]:
        """返回所有属性的 JSON 可序列化字典。"""
        raise NotImplementedError("Project.json")


__all__ = ["Project"]
