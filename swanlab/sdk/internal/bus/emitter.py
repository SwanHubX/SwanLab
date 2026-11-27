"""
@author: cunyue
@file: emitter.py
@time: 2026/3/13
@description: SwanLab 内部事件发射器，封装队列写入，仅供 SDK 内部组件使用
"""

import queue
from typing import Callable, Optional, Protocol

from swanlab.sdk.internal.pkg import console, throttle

from .events import EventPayload

__all__ = ["RunEmitter", "RunQueue", "EmitterProtocol"]

_DROP_WARNING_INTERVAL: float = 10.0
_QUEUE_PUT_POLL_INTERVAL: float = 0.5


RunQueue = queue.Queue[EventPayload]


class EmitterProtocol(Protocol):
    def emit(self, event: EventPayload) -> None: ...
    @property
    def queue(self) -> RunQueue: ...


class RunEmitter(EmitterProtocol):
    """
    内部事件发射器。

    持有事件队列，是唯一的写入入口。
    Run 公开 API 经过验证后调用 emit()；
    硬件监控、元数据采集等系统组件在构造时注入此对象，直接调用 emit()，
    不经过任何公开 API 或全局状态。
    """

    def __init__(self, maxsize: int = 100_000):
        self._queue: RunQueue = RunQueue(maxsize=maxsize)
        self._consumer_guard: Optional[Callable[[], bool]] = None
        self._drop_throttle = throttle.WarningThrottle(_DROP_WARNING_INTERVAL)

    @property
    def queue(self) -> RunQueue:
        """只读暴露给 BackgroundConsumer，外部不应直接操作队列"""
        return self._queue

    def set_consumer_guard(self, guard: Callable[[], bool]) -> None:
        """绑定消费者存活检查；由组件工厂注入。"""
        self._consumer_guard = guard

    def emit(self, event: EventPayload) -> None:
        """消费者存活时阻塞等待队列空位；线程已退出时丢弃，避免卡死训练线程。"""
        guard = self._consumer_guard
        if guard is None:
            self._queue.put(event, block=True)
            return
        while True:
            if not guard():
                self._warn_dropped()
                return
            try:
                self._queue.put(event, timeout=_QUEUE_PUT_POLL_INTERVAL)
                return
            except queue.Full:
                continue

    def _warn_dropped(self) -> None:
        count = self._drop_throttle.record()
        if count > 0:
            with console.skip_capture():
                console.warning(f"Dropped {count} event(s): background consumer is not running")
