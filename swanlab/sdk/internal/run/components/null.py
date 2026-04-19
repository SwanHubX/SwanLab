"""
@author: cunyue
@file: null.py
@time: 2026/4/19 19:02
@description: 空运行时组件，用于禁用运行时
"""

from queue import Queue

from swanlab.sdk.internal.bus.emitter import EmitterProtocol, RunQueue
from swanlab.sdk.internal.bus.events import EventPayload

from .consumer import ConsumerProtocol


class NullEmitter(EmitterProtocol):
    """空发射器，emit 为 no-op"""

    def __init__(self):
        self._queue = Queue()

    def emit(self, event: EventPayload) -> None:
        pass

    @property
    def queue(self) -> RunQueue:
        return self._queue


class NullConsumer(ConsumerProtocol):
    """空消费者，start/join 为 no-op"""

    def start(self) -> None:
        pass

    def stop(self) -> None:
        pass

    def join(self) -> None:
        pass
