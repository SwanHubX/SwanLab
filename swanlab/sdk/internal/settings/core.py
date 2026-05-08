"""
@author: cunyue
@file: core.py
@time: 2026/5/8
@description: SwanLab Core 配置
"""

import os

from pydantic import BaseModel, ConfigDict, Field


def section_rule_index_factory() -> int:
    # 向下兼容旧版本环境变量
    return int(os.environ.get("SWANLAB_CORE_SECTION_RULE_INDEX", os.environ.get("SWANLAB_SECTION_RULE_IDX", "0")))


class CoreSettings(BaseModel):
    section_rule_index: int = Field(default_factory=section_rule_index_factory)
    """
    用于将 metric key 拆分为 section name 的斜杠索引。
    0 表示第一个斜杠，1 表示第二个斜杠，-1 表示最后一个斜杠。
    超范围时会按 ``idx % slash_count`` 环绕。
    每次 run 的拆分位置会在 settings 加载时固定。
    """

    model_config = ConfigDict(frozen=True)
