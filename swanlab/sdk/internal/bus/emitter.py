"""
@author: cunyue
@file: emitter.py
@time: 2026/3/13
@description: SwanLab 内部事件发射器，封装队列写入，仅供 SDK 内部组件使用
"""

import queue

from .events import EventPayload

__all__ = ["RunEmitter"]


class RunEmitter:
    """
    内部事件发射器。

    持有事件队列，是唯一的写入入口。
    SwanLabRun 公开 API 经过验证后调用 emit()；
    硬件监控、元数据采集等系统组件在构造时注入此对象，直接调用 emit()，
    不经过任何公开 API 或全局状态。
    """

    def __init__(self, maxsize: int = 100_000):
        self._queue: queue.Queue[EventPayload] = queue.Queue(maxsize=maxsize)

    @property
    def queue(self) -> "queue.Queue[EventPayload]":
        """只读暴露给 BackgroundConsumer，外部不应直接操作队列"""
        return self._queue

    def emit(self, event: EventPayload) -> None:
        """将事件推入队列，背压时阻塞调用方"""
        self._queue.put(event, block=True)
