"""
@author: cunyue
@file: __init__.py
@time: 2026/4/11 18:45
@description: 异常处理模块，确保所有异常都被捕获，提供两个API
1. safe: 以装饰器的形式装饰一个函数，在函数内部捕获所有异常并返回None
2. safe_block: 以 with 的形式使用这个装饰器，在 with 语句块中捕获异常

在设计上这两个API是对console.trace的封装，console.trace会追踪当前调用栈，在终端打印简短信息，在日志文件中打印详细信息，并将详细信息写入日志文件
"""

import sys
from contextlib import contextmanager
from functools import wraps
from typing import Callable, Generator, Literal, Optional, Type, TypeVar

from swanlab.sdk.internal.pkg import console

__all__ = ["safe", "safe_block"]

# Python 3.10+ 直接用: from typing import ParamSpec
if sys.version_info >= (3, 10):
    from typing import ParamSpec
else:
    from typing_extensions import ParamSpec

# 定义参数规范 P 和返回值类型 R
P = ParamSpec("P")
R = TypeVar("R")


def safe(
    *exceptions: Type[BaseException],
    level: Literal["debug", "error"] = "error",
    message: Optional[str],
    write: bool = True,
    on_error: Optional[Callable[[BaseException], None]] = None,
):
    """
    捕获指定异常并返回 None
    使用 ParamSpec 和 TypeVar 动态修改类型签名，使 IDE 能正确推断出 Optional[R]。

    Args:
        *exceptions: 要捕获的异常类型，默认捕获所有 Exception
        level: 错误日志的级别，默认为error，可选debug，遵循console本身的level等级机制，如果为debug且不在debug模式，message和write都无效
        message: 错误信息，必传，如果传入None，则不打印错误信息
        write: 是否将错误信息写入日志，默认为 True，如果message为None，无论write是否为True，都不写入日志
        on_error: 异常发生时执行的回调函数，默认为 None
    """
    catch_types = exceptions if exceptions else (Exception,)

    # 这里的关键：明确告诉 IDE，接收 Callable[P, R]，返回 Callable[P, Optional[R]]
    def decorator(func: Callable[P, R]) -> Callable[P, Optional[R]]:

        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> Optional[R]:
            try:
                return func(*args, **kwargs)
            except catch_types as e:
                if message is not None:
                    console.trace(message, write_to_file=write, level_name=level)
                if on_error is not None:
                    on_error(e)

        return wrapper

    return decorator


@contextmanager
def safe_block(
    *exceptions: Type[BaseException],
    level: Literal["debug", "error"] = "error",
    message: Optional[str],
    write: bool = True,
    on_error: Optional[Callable[[BaseException], None]] = None,
) -> Generator[None, None, None]:
    """
    以 with 语句的形式捕获指定异常，异常发生时打印 trace 日志后静默退出块。

    Args:
        *exceptions: 要捕获的异常类型，默认捕获所有 Exception
        level: 错误日志的级别，默认为error，可选debug，遵循console本身的level等级机制，如果为debug且不在debug模式，message和write都无效
        message: 错误信息，必传，如果传入None，则不打印错误信息
        write: 是否将错误信息写入日志，默认为 True
        on_error: 异常发生时执行的回调函数，默认为 None

    Example:
        with safe_block(message="failed to do something"):
            risky_operation()
    """
    catch_types = exceptions if exceptions else (Exception,)
    try:
        yield
    except catch_types as e:
        if message is not None:
            console.trace(message, write_to_file=write, level_name=level)
        if on_error is not None:
            on_error(e)
