"""
@author: cunyue
@file: __init__.py
@time: 2026/3/7 14:18
@description: SwanLab SDK 辅助函数
"""

import sys
from functools import wraps
from typing import Callable, Optional, Type, TypeVar

from . import rich

__all__ = ["catch_and_return_none", "rich"]

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
