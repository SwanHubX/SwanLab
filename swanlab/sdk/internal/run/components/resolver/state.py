"""
@author: caddiesnew
@file: state.py
@time: 2026/8/10
@description: define_metric resolver 状态数据结构。
"""

import sys
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

__all__ = [
    "EffectiveDefinition",
    "DefinitionPatch",
    "RuleState",
    "ConcreteState",
    "make_effective",
]

# py3.9 兼容：slots 仅 3.10+ 支持
# 3.9 退化为普通 dict-based dataclass
_SLOTS_KWARGS = {"slots": True} if sys.version_info >= (3, 10) else {}


@dataclass(frozen=True, **_SLOTS_KWARGS)
class EffectiveDefinition:
    """merge/replace 后的完整 effective 快照。"""

    x_axis: str  # X 轴 metric key；"_step" 表示系统默认
    section_name: Optional[str]  # 图表分组名称；None 表示使用默认 section
    hidden: bool  # 是否隐藏该 metric 对应的图表
    step_sync: bool = True  # 该 Y 是否在 X 分次 log 时触发 custom X 注入；默认开启


# EffectiveDefinition 是 frozen+hashable，相同 (x_axis, section_name, hidden, step_sync)
# 复用同一实例。intern cache 为进程级，distinct 组合数量有界，无需清理。
_EFF_INTERN: Dict[Tuple[str, Optional[str], bool, bool], EffectiveDefinition] = {}


def make_effective(
    x_axis: str, section_name: Optional[str], hidden: bool, step_sync: bool = True
) -> EffectiveDefinition:
    """构造并 intern 一个 EffectiveDefinition。"""
    key = (x_axis, section_name, hidden, step_sync)
    eff = _EFF_INTERN.get(key)
    if eff is None:
        eff = EffectiveDefinition(x_axis=x_axis, section_name=section_name, hidden=hidden, step_sync=step_sync)
        _EFF_INTERN[key] = eff
    return eff


@dataclass(**_SLOTS_KWARGS)
class DefinitionPatch:
    """presence-aware 定义补丁。

    None 表示 "未提供"，merge 时保留旧值，replace 时重置为默认。
    """

    x_axis: Optional[str] = None  # X 轴补丁；None 表示调用方未提供
    section_name: Optional[str] = None  # section 补丁；None 表示调用方未提供
    hidden: Optional[bool] = None  # 隐藏状态补丁；None 表示调用方未提供
    step_sync: Optional[bool] = None  # 注入开关补丁；None 表示调用方未提供，True/False 为显式值


@dataclass(**_SLOTS_KWARGS)
class RuleState:
    """一条 define_metric rule（exact 或 glob）的状态。"""

    identity: str  # exact key 或 glob pattern，用于标识 rule
    is_glob: bool  # 是否为末尾带 * 的 prefix glob rule
    effective: EffectiveDefinition  # rule 经 merge/replace 后的完整定义快照
    patch: DefinitionPatch  # 本次 define_metric 调用提供的 presence-aware 补丁
    overwrite: bool  # True 表示从默认值 replace，False 表示在已有状态上 merge


@dataclass(**_SLOTS_KWARGS)
class ConcreteState:
    """一个已物化的 concrete key 的状态。"""

    key: str  # log 中实际出现的 metric key
    metric_class: str  # resolver identity 的类别：SCALAR 或 MEDIA
    effective: EffectiveDefinition  # 首次物化或 exact 更新后固定的定义快照
    source: str  # 定义来源：exact、glob 或 automatic
    emitted_effective: Optional[EffectiveDefinition] = None  # 最近一次已生成 ColumnRecord 的定义快照
    column_type: Optional[int] = None  # 已识别的 protobuf ColumnType；未见数据时为 None
