"""
@author: caddiesnew
@file: common.py
@time: 2026/4/20
@description: 公共查询 API 通用类型定义
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, TypedDict

# 启用/停用
ApiStatusLiteral = Literal["ENABLED", "DISABLED"]

# 侧边列类型
# STABLE: Experiment 的固有字段，如 state, name, labels, colors 等
# CONFIG: 动态生成的实验配置字段，如 learning_rate, batch_size 等
# SCALAR: 动态生成的标量字段，一般用于标量图展示，如 train/loss 等
ApiSidebarLiteral = Literal["SCALAR", "CONFIG", "STABLE"]

# 实验类型： 运行中/总览
ApiExperimentTypeLiteral = Literal["CHAPTER", "SUMMARY"]

# 实验状态类型
ApiRunStateLiteral = Literal["RUNNING", "FINISHED", "CRASHED", "ABORTED", "OFFLINE"]

# 可见性类型
ApiVisibilityLiteral = Literal["PUBLIC", "PRIVATE"]

# 工作空间类型
ApiWorkspaceLiteral = Literal["TEAM", "PERSON"]

# 工作空间成员类型
ApiRoleLiteral = Literal["VISITOR", "VIEWER", "MEMBER", "OWNER"]

# 列种类
ApiColumnClassLiteral = Literal["CUSTOM", "SYSTEM"]
# 列数据类型
ApiColumnDataTypeLiteral = Literal[
    "FLOAT",
    "BOOLEAN",
    "STRING"
    # media 类型
    "IMAGE",
    "AUDIO",
    "VIDEO",
    # 3D点云 (json)
    "OBJECT3D",
    #   生物化学分子
    "MOLECULE",
    # (js/ ts 文件)
    "ECHARTS",
    # 表格类型
    "TABLE",
    "TEXT",
]

# Self-Hosted 身份类型
ApiIdentityLiteral = Literal["root", "user"]

# License 许可证类型
ApiLicensePlanLiteral = Literal["free", "commercial"]


# 后端允许的每页条数
_VALID_PAGE_SIZES = (10, 12, 15, 20, 24, 27, 50, 100)


@dataclass(frozen=True)
class PaginatedQuery:
    """
    通用分页查询参数，与后端 pagination_query 对齐。

    page: 当前页码，≥1
    size: 每页条数，必须为后端允许值之一
    search: 搜索关键词
    sort: 排序字段
    all: 是否拉取全部分页（客户端侧自动翻页）
    """

    page: int = 1
    size: int = 20
    search: Optional[str] = None
    sort: Optional[str] = None
    all: bool = False

    def __post_init__(self) -> None:
        if self.page < 1:
            raise ValueError(f"page must be >= 1, got {self.page}")
        if self.size not in _VALID_PAGE_SIZES:
            raise ValueError(f"size must be one of {_VALID_PAGE_SIZES}, got {self.size}")

    def to_params(self, **extra: Optional[Any]) -> Dict[str, Any]:
        """转换为查询参数字典，自动过滤 None 值。"""
        params: Dict[str, Any] = {"page": self.page, "size": self.size}
        if self.search is not None:
            params["search"] = self.search
        if self.sort is not None:
            params["sort"] = self.sort
        params.update({k: v for k, v in extra.items() if v is not None})
        return params


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


__all__ = ["ApiPaginationType", "ApiResponseType", "PaginatedQuery"]
