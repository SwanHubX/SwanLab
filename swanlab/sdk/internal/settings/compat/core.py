"""Core 域遗留环境变量兼容层。"""

from .env import getenv

__all__ = ["section_rule_index_factory"]


def section_rule_index_factory() -> int:
    # 向下兼容旧版本环境变量
    return int(getenv("SWANLAB_SECTION_RULE_IDX", "0"))
