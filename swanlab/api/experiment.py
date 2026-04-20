"""
@author: Nexisato
@file: experiment.py
@time: 2026/4/20
@description: Experiment 实体类 — 单个实验的查询与操作
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional

from swanlab.sdk.typings.core_python.api.experiment import RunType

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class Experiment:
    """
    表示一个 SwanLab 实验。

    通过独立 Client 实例与云端交互，支持查询实验详情、获取指标、删除等操作。
    """

    def __init__(self, client: "Client", *, path: str, web_host: str) -> None:
        self._client: "Client" = client
        self._path: str = path  # 'username/project/run_id'
        self._web_host: str = web_host
        self._data: Optional[RunType] = None

    def _fetch(self) -> RunType:
        """从云端拉取实验详情，缓存到 _data。"""
        raise NotImplementedError("Experiment._fetch")

    # ------------------------------------------------------------------
    #  属性占位
    # ------------------------------------------------------------------

    @property
    def id(self) -> str:
        """实验 CUID（唯一标识）。"""
        raise NotImplementedError

    @property
    def name(self) -> str:
        """实验名称。"""
        raise NotImplementedError

    @property
    def state(self) -> str:
        """实验状态（RUNNING / FINISHED / CRASHED / ABORTED）。"""
        raise NotImplementedError

    @property
    def url(self) -> str:
        """实验在 Web 面板中的访问 URL。"""
        raise NotImplementedError

    # ------------------------------------------------------------------
    #  操作占位
    # ------------------------------------------------------------------

    def metrics(self, keys: Optional[List[str]] = None) -> Any:
        """获取实验指标数据。"""
        raise NotImplementedError("Experiment.metrics")

    def delete(self) -> None:
        """删除此实验。"""
        raise NotImplementedError("Experiment.delete")

    def json(self) -> Dict[str, Any]:
        """返回所有属性的 JSON 可序列化字典。"""
        raise NotImplementedError("Experiment.json")


__all__ = ["Experiment"]
