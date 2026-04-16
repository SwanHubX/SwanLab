"""
@author: cunyue
@file: core.py
@time: 2026/3/13
@description: SwanLab Core 接口协议，负责对产出的Record做持久化处理和后端交互
"""

import threading
from abc import ABC, abstractmethod
from enum import Enum
from typing import Generic, List, Optional, TypeVar, cast

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.context import RunContext

__all__ = ["CoreEnum", "CoreProtocol", "DeliverHandle"]


class CoreEnum(str, Enum):
    CORE_PYTHON = "CorePython"
    CORE = "Core"


T = TypeVar("T")

_UNSET = object()
"""
哨兵值，解决泛型类型推导问题
"""


class DeliverHandle(Generic[T]):
    """deliver 操作的确认句柄。

    调用方通过 wait() 阻塞等待云端确认；
    实现方通过 resolve() / reject() 完成确认。
    """

    def __init__(self):
        self._event = threading.Event()
        self._result: object = _UNSET
        self._exception: Optional[BaseException] = None

    def wait(self, timeout: Optional[float] = None) -> T:
        """阻塞等待确认，返回结果。超时抛 TimeoutError，处理失败抛原始异常。"""
        if not self._event.wait(timeout=timeout):
            raise TimeoutError("Deliver timed out, please check your network connection")
        if self._exception is not None:
            raise self._exception
        assert self._result is not _UNSET  # resolve 必须在 _event.set() 之前调用
        return cast(T, self._result)

    def resolve(self, result: T = None) -> None:
        """确认成功，设置结果并唤醒等待方。"""
        self._result = result
        self._event.set()

    def reject(self, exc: BaseException) -> None:
        """确认失败，设置异常并唤醒等待方。"""
        self._exception = exc
        self._event.set()


class CoreProtocol(ABC):
    """
    SwanLab Core 协议
    协议层作为一个“垫片”，抹平前端与后端实现的差异
    """

    def __init__(self, ctx: RunContext):
        self._ctx = ctx

    @abstractmethod
    def publish(self, records: List[Record]) -> None:
        """即发即忘：持久化 + 推上传队列，不等待确认。用于高频数据。"""
        ...

    @abstractmethod
    def deliver(self, record: Record) -> DeliverHandle:
        """
        等待确认：返回 DeliverHandle，在记录持久化+上传完成后 resolve。
        调用方可通过 handle.wait(timeout=30) 阻塞等待，设计上仅用于 run_start、finish 等生命周期事件
        """
        ...

    @abstractmethod
    def fork(self) -> "CoreProtocol":
        """
        创建一个新的 CoreProtocol，此函数用于处理多进程fork的情况
        """
        ...
