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
    "STRING",
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

# 列数据非 media 类型，方便过滤
ApiColumnScalarTypeLiteral = Literal["FLOAT", "BOOLEAN", "STRING"]

# Self-Hosted 身份类型
ApiIdentityLiteral = Literal["root", "user"]

# License 许可证类型
ApiLicensePlanLiteral = Literal["free", "commercial"]

# 指标类型（log 不属于 column-backed metrics，使用独立查询方法）
ApiMetricColumnTypeLiteral = Literal["SCALAR", "MEDIA"]

# 指标扩展类型（包含 LOG，用于内部 Metric 调度）
ApiMetricAllTypeLiteral = Literal["SCALAR", "MEDIA", "LOG"]

# 指标日志级别
ApiMetricLogLevelLiteral = Literal["DEBUG", "INFO", "WARN", "ERROR"]

# X 轴类型
ApiMetricXAxisLiteral = Literal["step", "time", "relative_time"]


# ---------------------------------------------------------------------------
# STABLE 字段 key 枚举
# 对应 experiment 表的直接字段或嵌套字段，来自 sidebar.js stableFieldSelect。
# ---------------------------------------------------------------------------
ApiFilterStableKeyLiteral = Literal[
    # 实验状态 RUNNING / FINISHED / CRASHED / ABORTED
    "state",
    # 实验名称
    "name",
    # 实验描述
    "description",
    # 是否可见
    "show",
    # TODO: 已知被设置为 pin 的 experiment 在任何情况下会强制返回
    "pin",
    # 是否为基线
    "baseline",
    # 颜色
    "colors",
    # 实验分组名
    "cluster",
    # 分布式任务类型
    "job",
    # 创建时间
    "createdAt",
    # 更新时间
    "updatedAt",
    # 完成时间
    "finishedAt",
    # 收藏时间
    "pinnedAt",
    # 标签名数组
    "labels",
]

# ---------------------------------------------------------------------------
# 过滤操作符
# ---------------------------------------------------------------------------
# EQ        : 等于
# NEQ       : 不等于
# GTE       : 大于等于（数值 / 日期 / 字符串）
# LTE       : 小于等于（数值 / 日期 / 字符串）
# IN        : 在给定值列表中
# NOT IN    : 不在给定值列表中
# CONTAIN   : 模糊包含
#
# 注意：
# - 数组类型（如 labels）仅支持 EQ / NEQ / IN / NOT IN / CONTAIN
# - 日期类型 GTE/LTE 用 Date 对象比较；其余用 ISO 字符串比较
# - 数值类型优先数值比较，失败回退字符串比较
# ---------------------------------------------------------------------------
ApiFilterOpLiteral = Literal["EQ", "NEQ", "GTE", "LTE", "IN", "NOT IN", "CONTAIN"]

# ---------------------------------------------------------------------------
# 排序方向
# ---------------------------------------------------------------------------
ApiSortOrderLiteral = Literal["ASC", "DESC"]

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
