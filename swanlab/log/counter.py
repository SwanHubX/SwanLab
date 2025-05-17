"""
@author: cunyue
@file: counter.py
@time: 2025/5/17 17:23
@description: 线程安全计数器，用于记录当前epoch
"""

import threading
from contextlib import ContextDecorator


class AtomicCounter(ContextDecorator):
    """A thread-safe counter that allows free increments within `with` context."""

    def __init__(self, initial_value: int = 0) -> None:
        if initial_value < 0:
            raise ValueError("Counter initial value cannot be negative.")
        self._value = initial_value
        self._lock = threading.Lock()
        self._context_entered = False  # 跟踪是否在上下文管理中

    @property
    def value(self) -> int:
        """Get the current value (thread-safe)."""
        with self._lock:
            return self._value

    def increment(self) -> int:
        """Increment the counter by 1 and return the new value."""
        if not self._context_entered:
            # 不在上下文中时报错
            raise RuntimeError("Counter increment is only allowed within the context manager.")
        else:
            # 在上下文中时直接操作（锁已在__enter__中获取）
            self._value += 1
            return self._value

    def __enter__(self) -> 'AtomicCounter':
        """Enter the context: acquire the lock for batch operations."""
        self._lock.acquire()
        self._context_entered = True
        return self  # 返回self以支持链式调用

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Exit the context: release the lock."""
        self._context_entered = False
        self._lock.release()

    def __str__(self) -> str:
        return f"AtomicCounter(value={self.value})"
