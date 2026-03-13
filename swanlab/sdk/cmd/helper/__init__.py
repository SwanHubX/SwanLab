"""
@author: cunyue
@file: __init__.py
@time: 2026/3/13 21:37
@description: SwanLab SDK CMD 辅助函数
"""

import threading
from functools import wraps

_CMD_LOCK = threading.Lock()


def with_cmd_lock(func):
    """
    全局锁装饰器。
    注意：此锁为不可重入锁 (Lock)，如果两个被此装饰器装饰的 API 相互调用，必定发生死锁
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        with _CMD_LOCK:
            return func(*args, **kwargs)

    return wrapper
