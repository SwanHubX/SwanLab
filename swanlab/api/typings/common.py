"""
@author: caddiesnew
@file: common.py
@time: 2026/4/20
@description: 公共查询 API 通用类型定义
"""

from typing import Any, List, TypedDict


class ApiLabelType(TypedDict):
    name: str


class ApiPaginationType(TypedDict):
    list: List
    size: int
    pages: int
    total: int


class ApiResponseType(TypedDict):
    """API 响应的统一封装类型。"""

    ok: bool
    errmsg: str
    data: Any


class ApiResponse:
    """
    API 响应的统一封装，保证任何异常都不会导致程序 crash。

    - ok=True  时 data 持有正常返回值
    - ok=False 时 data 为 None，errmsg 描述失败原因
    """

    __slots__ = ("ok", "errmsg", "data")

    def __init__(self, *, ok: bool, errmsg: str = "", data: Any = None) -> None:
        self.ok = ok
        self.errmsg = errmsg
        self.data = data

    def to_dict(self) -> ApiResponseType:
        return {"ok": self.ok, "errmsg": self.errmsg, "data": self.data}

    def to_json_dict(self) -> dict:
        """返回 JSON 可序列化的字典，自动将实体 data 转为 dict。"""
        data = self.data
        if data is not None and hasattr(data, "to_dict"):
            data = data.to_dict()
        return {"ok": self.ok, "errmsg": self.errmsg, "data": data}

    def __repr__(self) -> str:
        if self.ok:
            return f"ApiResponse(ok=True, data={self.data!r})"
        return f"ApiResponse(ok=False, errmsg={self.errmsg!r})"
