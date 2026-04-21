"""
@author: caddiesnew
@file: common.py
@time: 2026/4/20
@description: 公共查询 API 通用类型定义
"""

from typing import Any, Dict, List, Literal, TypedDict

# 实验状态类型
ApiRunStateType = Literal["RUNNING", "FINISHED", "CRASHED", "ABORTED", "OFFLINE"]

# 可见性类型
ApiVisibilityType = Literal["PUBLIC", "PRIVATE"]

# 工作空间类型
ApiWorkspaceType = Literal["TEAM", "PERSON"]

# 工作空间成员类型
ApiRoleType = Literal["VISITOR", "VIEWER", "MEMBER", "OWNER"]

# Self-Hosted 身份类型
ApiIdentityType = Literal["root", "user"]

# License 许可证类型
ApiLicensePlanType = Literal["free", "commercial"]


class ApiLabelType(TypedDict):
    name: str


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

    def to_dict(self) -> Dict[str, Any]:
        return {"ok": self.ok, "errmsg": self.errmsg, "data": self.data}

    def to_json_dict(self) -> Dict[str, Any]:
        """返回 JSON 可序列化的字典，自动将实体 data 转为 dict。"""
        data = self.data
        errors: list[str] = []
        if not self.ok and self.errmsg:
            errors.append(self.errmsg)
        if data is not None and hasattr(data, "to_dict"):
            data = data.to_dict()
            # 收集实体内部子请求的错误
            if hasattr(data, "__getitem__"):
                # to_dict 返回的 dict 不带 _errors，需要从实体取
                pass
            entity_errors = getattr(self.data, "_errors", [])
            errors.extend(entity_errors)
        ok = self.ok and not errors
        return {"ok": ok, "errmsg": "; ".join(errors) if errors else "", "data": data}

    def __repr__(self) -> str:
        if self.ok:
            return f"ApiResponse(ok=True, data={self.data!r})"
        return f"ApiResponse(ok=False, errmsg={self.errmsg!r})"


__all__ = ["ApiLabelType", "ApiPaginationType", "ApiResponseType"]
