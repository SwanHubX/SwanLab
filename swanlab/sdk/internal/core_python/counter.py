"""
@author: cunyue
@file: counter.py
@time: 2026/4/21 19:04
@description: 线程安全计数器
"""

from __future__ import annotations

import threading


class Counter:
    """线程安全计数器。"""

    def __init__(self, initial: int = 0) -> None:
        self._value = initial
        self._lock = threading.Lock()

    @property
    def value(self) -> int:
        """获取当前计数值。"""
        with self._lock:
            return self._value

    def inc(self, delta: int = 1) -> int:
        """计数增加，返回增加后的值。"""
        with self._lock:
            self._value += delta
            return self._value

    def dec(self, delta: int = 1) -> int:
        """计数减少，返回减少后的值。"""
        with self._lock:
            self._value -= delta
            return self._value

    def reset(self, value: int = 0) -> None:
        """重置计数器。"""
        with self._lock:
            self._value = value

    def get_and_reset(self, value: int = 0) -> int:
        """获取当前值并重置。"""
        with self._lock:
            old = self._value
            self._value = value
            return old

    def __int__(self) -> int:
        return self.value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(value={self.value})"
