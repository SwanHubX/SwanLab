"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:18
@description: SwanLab SDK 辅助函数
"""

import sys
from functools import wraps
from typing import Callable, Optional, Type, TypeVar

from . import env, rich

__all__ = ["catch_and_return_none", "rich", "env", "strip_none"]

# Python 3.10+ 直接用: from typing import ParamSpec
if sys.version_info >= (3, 10):
    from typing import ParamSpec
else:
    from typing_extensions import ParamSpec

# 定义参数规范 P 和返回值类型 R
P = ParamSpec("P")
R = TypeVar("R")


def catch_and_return_none(*exceptions: Type[BaseException]):
    """
    捕获指定异常并返回 None。
    使用 ParamSpec 和 TypeVar 动态修改类型签名，使 IDE 能正确推断出 Optional[R]。
    """
    catch_types = exceptions if exceptions else (Exception,)

    # 这里的关键：明确告诉 IDE，接收 Callable[P, R]，返回 Callable[P, Optional[R]]
    def decorator(func: Callable[P, R]) -> Callable[P, Optional[R]]:

        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> Optional[R]:
            try:
                return func(*args, **kwargs)
            except catch_types:
                return None

        return wrapper

    return decorator


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
