"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13 21:37
@description: SwanLab SDK CMD 辅助函数
"""

import os
import threading
from functools import wraps
from typing import Callable

from swanlab.sdk.internal.run import has_run

_CMD_LOCK = threading.Lock()
_CMD_PID = os.getpid()


def with_cmd_lock(func):
    """
    全局锁装饰器。
    注意：此锁为不可重入锁 (Lock)，如果两个被此装饰器装饰的 API 相互调用，必定发生死锁
    fork 后子进程中此锁可能处于已持有状态，检测到 fork 时替换为新锁
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        global _CMD_LOCK, _CMD_PID
        if os.getpid() != _CMD_PID:
            _CMD_LOCK = threading.Lock()
            _CMD_PID = os.getpid()
        with _CMD_LOCK:
            return func(*args, **kwargs)

    return wrapper


def with_run(cmd: str):
    """
    装饰器：要求必须有 run 在运行，否则抛出 RuntimeError

    :param cmd: 命令名称
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            if not has_run():
                raise RuntimeError(f"`swanlab.{cmd}` requires an active Run, call `swanlab.init()` first.")
            return func(*args, **kwargs)

        return wrapper

    return decorator


def without_run(cmd: str):
    """
    装饰器：要求必须没有 run 在运行，否则抛出 RuntimeError
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            if has_run():
                raise RuntimeError(f"`swanlab.{cmd}` requires no active Run, call `swanlab.finish()` first.")
            return func(*args, **kwargs)

        return wrapper

    return decorator
