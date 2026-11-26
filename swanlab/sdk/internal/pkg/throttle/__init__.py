"""
@author: Nexisato
@file: __init__.py
@time: 2026/11/26
@description: 基于单调时钟的轻量级告警节流工具，供事件队列、消费者与状态监听复用。
"""

import threading
import time
from typing import Callable, Optional

__all__ = ["WarningThrottle"]


class WarningThrottle:
    """线程安全的单调时钟告警节流器。"""

    def __init__(self, interval: float, *, clock: Callable[[], float] = time.monotonic) -> None:
        self._interval = interval
        self._clock = clock
        self._lock = threading.Lock()
        self._last_at: Optional[float] = None
        self._pending: int = 0

    def should_emit(self) -> bool:
        """若当前处于可打印窗口则推进窗口并返回 True，否则返回 False。"""
        now = self._clock()
        with self._lock:
            if self._last_at is None or now - self._last_at >= self._interval:
                self._last_at = now
                return True
            return False

    def record(self, count: int = 1) -> int:
        """累计待报告次数；到达打印窗口时清零并返回累计总数，否则返回 0。"""
        now = self._clock()
        with self._lock:
            self._pending += count
            if self._last_at is None or now - self._last_at >= self._interval:
                total = self._pending
                self._pending = 0
                self._last_at = now
                return total
            return 0

    def reset(self) -> None:
        """重置节流窗口与累计计数，下一次触发立即允许打印。"""
        with self._lock:
            self._last_at = None
            self._pending = 0
