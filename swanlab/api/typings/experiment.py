"""
@author: caddiesnew
@file: experiment.py
@time: 2026/4/20
@description: 公共查询 API 实验类型定义

POST /runs/shows 接口支持三个维度的筛选和组织：过滤 (filters)、分组 (groups)、排序 (sorts)。
每项都有一个 type 字段，取值取决于数据来源：
  - STABLE: 实验表固有字段（固定枚举）
  - CONFIG: experimentProfile.config 中用户定义的超参（动态 key，如 learning_rate）
  - SCALAR: 训练过程中记录的标量指标最新值（动态 key，如 train/loss）
"""

from typing import Any, Dict, List, Literal, Optional, TypedDict

from .common import ApiExperimentTypeLiteral, ApiRunStateLiteral, ApiSidebarLiteral
from .user import ApiUserType

# ---------------------------------------------------------------------------
# STABLE 字段 key 枚举
# 对应 experiment 表的直接字段或嵌套字段，来自 sidebar.js stableFieldSelect。
# ---------------------------------------------------------------------------
ApiStableKeyLiteral = Literal[
    # 实验状态 RUNNING / FINISHED / CRASHED / ABORTED
    "state",
    # 实验名称
    "name",
    # 实验描述
    "description",
    # 是否可见
    "show",
    # TODO: experiment 被设置为 pin 时强制返回
    # "pin",
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


# ---------------------------------------------------------------------------
# filter / group / sort item
# POST /runs/shows 请求体中 filters / groups / sorts 数组的元素类型。
# 用户传入时不需要 active 字段，由 SDK 内部自动补充 active: True。
# ---------------------------------------------------------------------------
class ApiFilterItem(TypedDict):
    """POST /runs/shows 请求体中的过滤项。

    多个 filter 之间为 AND 关系。
    收藏的实验（pin: true）永远不会被过滤掉。
    """

    key: str  # STABLE 时为 ApiStableKeyLiteral 枚举；CONFIG/SCALAR 时为动态字段名
    type: ApiSidebarLiteral  # STABLE | CONFIG | SCALAR
    op: ApiFilterOpLiteral  # 过滤操作符
    value: List[str]  # 过滤值列表（空值统一视为空字符串 ""）


class ApiGroupItem(TypedDict):
    """POST /runs/shows 请求体中的分组项。

    多个 group 形成多层嵌套（外层 group 为第一层）。
    数组类型值（如 labels）会排序后用 ", " 连接成字符串作为分组 key。
    """

    key: str  # STABLE 时为 ApiStableKeyLiteral 枚举；CONFIG/SCALAR 时为动态字段名
    type: ApiSidebarLiteral  # STABLE | CONFIG | SCALAR


class ApiSortItem(TypedDict):
    """POST /runs/shows 请求体中的排序项。

    排序与分组联动：有 order 的字段按方向排序后平铺，无 order 的保留嵌套结构。
    后端自动追加兜底排序：pin DESC > pinnedAt DESC > createdAt DESC。
    """

    key: str  # STABLE 时为 ApiStableKeyLiteral 枚举；CONFIG/SCALAR 时为动态字段名
    type: ApiSidebarLiteral  # STABLE | CONFIG | SCALAR
    order: ApiSortOrderLiteral  # ASC | DESC


# ---------------------------------------------------------------------------
# 实验实体
# ---------------------------------------------------------------------------
class ApiExperimentLabelType(TypedDict):
    name: str


class ApiExperimentProfileType(TypedDict):
    config: Dict[str, Any]
    metadata: Dict[str, Any]
    requirements: str
    conda: str


class ApiExperimentType(TypedDict, total=False):
    project_id: str
    cuid: str
    name: str
    type: ApiExperimentTypeLiteral
    description: str
    labels: List[ApiExperimentLabelType]
    profile: ApiExperimentProfileType
    show: bool
    state: ApiRunStateLiteral
    cluster: str
    job: str
    user: ApiUserType
