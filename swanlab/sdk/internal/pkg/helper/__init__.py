"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:18
@description: SwanLab SDK 辅助函数
"""

from .env import DEBUG, is_interactive, is_jupyter
from .rich import with_loading_animation
from .version import get_swanlab_latest_version, get_swanlab_version

__all__ = [
    "DEBUG",
    "is_jupyter",
    "is_interactive",
    "strip_none",
    "with_loading_animation",
    "get_swanlab_version",
    "get_swanlab_latest_version",
]


def strip_none(data: dict, strip_empty_dict: bool = True, strip_empty_str: bool = False) -> dict:
    """
    递归剔除字典中的 None 值。

    :param data: 待处理的字典
    :param strip_empty_dict: 是否连同空字典一起剔除（包含剔除 None 后变为空的嵌套字典），默认为 True
    :param strip_empty_str: 是否连同空字符串一起剔除，仅当 strip_empty_dict 为 True 时有效，默认为 False
    """
    clean_data = {}
    for k, v in data.items():
        if isinstance(v, dict):
            # 递归调用时，记得把参数透传下去
            cleaned_v = strip_none(v, strip_empty_dict=strip_empty_dict, strip_empty_str=strip_empty_str)

            # 只有当嵌套字典里有值，或者我们明确声明“不剔除空字典”时，才保留这个 key
            if cleaned_v or not strip_empty_dict:
                clean_data[k] = cleaned_v
        elif v is not None and (not strip_empty_str or v != ""):
            clean_data[k] = v

    return clean_data
