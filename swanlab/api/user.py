"""
@author: Nexisato
@file: user.py
@time: 2026/4/20
@description: User 实体类 — 用户信息与 API Key 管理
"""

from typing import TYPE_CHECKING, Any, Dict, List

from swanlab.sdk.typings.core_python.api.user import ApiKeyType

if TYPE_CHECKING:
    from swanlab.sdk.internal.pkg.client import Client


class User:
    """
    表示一个 SwanLab 用户。

    通过独立 Client 实例与云端交互，支持查询用户信息、管理 API Key。
    """

    def __init__(self, client: "Client", *, username: str) -> None:
        self._client: "Client" = client
        self._username: str = username

    # ------------------------------------------------------------------
    #  属性占位
    # ------------------------------------------------------------------

    @property
    def username(self) -> str:
        """用户名。"""
        raise NotImplementedError

    @property
    def teams(self) -> List[str]:
        """用户所属的团队列表。"""
        raise NotImplementedError

    # ------------------------------------------------------------------
    #  操作占位
    # ------------------------------------------------------------------

    def api_keys(self) -> List[ApiKeyType]:
        """获取用户的 API Key 列表。"""
        raise NotImplementedError("User.api_keys")

    def json(self) -> Dict[str, Any]:
        """返回所有属性的 JSON 可序列化字典。"""
        raise NotImplementedError("User.json")


__all__ = ["User"]
