"""
@author: caddiesnew
@file: state.py
@time: 2026/8/10
@description: define_metric resolver 状态数据结构。
"""

import sys
from dataclasses import dataclass
from typing import Optional

__all__ = [
    "EffectiveDefinition",
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
    step_sync: bool = False  # 该 Y 是否在 X 分次 log 时触发 custom X 注入；默认值 = is_custom_x(x_axis)，显式指定优先


def make_effective(
    x_axis: str, section_name: Optional[str], hidden: bool, step_sync: bool = True
) -> EffectiveDefinition:
    """构造一个 EffectiveDefinition。"""
    return EffectiveDefinition(x_axis=x_axis, section_name=section_name, hidden=hidden, step_sync=step_sync)


@dataclass(**_SLOTS_KWARGS)
class ConcreteState:
    """一个已物化的 concrete key 的状态。"""

    key: str  # log 中实际出现的 metric key
    metric_class: str  # resolver identity 的类别：SCALAR 或 MEDIA
    effective: EffectiveDefinition  # 首次物化或 exact 更新后固定的定义快照
    source: str  # 定义来源：exact、glob；automatic 表示无 define 下被 log（仅钉住，不产出列）
    emitted: bool = False  # 是否已产出过 ColumnRecord（定义不回溯，每 (class, key) 只发一次）
