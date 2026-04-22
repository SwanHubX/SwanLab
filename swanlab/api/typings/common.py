"""
@author: caddiesnew
@file: common.py
@time: 2026/4/20
@description: 公共查询 API 通用类型定义
"""

from typing import Any, Dict, List, Literal, TypedDict

# 列类型
ApiColumnLiteral = Literal["SCALAR", "CONFIG", "STABLE"]

# 实验状态类型
ApiRunStateLiteral = Literal["RUNNING", "FINISHED", "CRASHED", "ABORTED", "OFFLINE"]

# 可见性类型
ApiVisibilityLiteral = Literal["PUBLIC", "PRIVATE"]

# 工作空间类型
ApiWorkspaceLiteral = Literal["TEAM", "PERSON"]

# 工作空间成员类型
ApiRoleLiteral = Literal["VISITOR", "VIEWER", "MEMBER", "OWNER"]

# Self-Hosted 身份类型
ApiIdentityLiteral = Literal["root", "user"]

# License 许可证类型
ApiLicensePlanLiteral = Literal["free", "commercial"]


class ApiPaginationType(TypedDict):
    list: List
    size: int
    pages: int
    total: int


class ApiResponseType:
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

    def json(self) -> Dict[str, Any]:
        """返回 JSON 可序列化的字典，自动将实体 data 转为 dict。"""
        data = self.data
        errors: list[str] = []
        if not self.ok and self.errmsg:
            errors.append(self.errmsg)
        if data is not None and hasattr(data, "json"):
            data = data.json()
            entity_errors = getattr(self.data, "_errors", [])
            errors.extend(entity_errors)
        ok = self.ok and not errors
        return {"ok": ok, "errmsg": "; ".join(errors) if errors else "", "data": data}

    def __repr__(self) -> str:
        if self.ok:
            return f"ApiResponse(ok=True, data={self.data!r})"
        return f"ApiResponse(ok=False, errmsg={self.errmsg!r})"


__all__ = ["ApiPaginationType", "ApiResponseType"]
