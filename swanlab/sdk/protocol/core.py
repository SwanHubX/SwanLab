"""
@author: cunyue
@file: core.py
@time: 2026/3/13
@description: SwanLab Core 接口协议，负责对产出的Record做持久化处理和后端交互
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import TYPE_CHECKING, List

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.run.v1.run_pb2 import (
    FinishRecord,
    FinishResponse,
    StartRecord,
    StartResponse,
)

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext

__all__ = ["CoreEnum", "CoreProtocol"]


class CoreEnum(str, Enum):
    CORE_PYTHON = "CorePython"
    CORE = "Core"


class CoreProtocol(ABC):
    """
    SwanLab Core 协议
    协议层作为一个“垫片”，抹平前端与后端实现的差异
    """

    def __init__(self, ctx: "RunContext"):
        self._ctx = ctx

    @abstractmethod
    def deliver_run_start(self, start_record: StartRecord) -> StartResponse:
        """
        交付运行开始事件，同步事件，等待确认
        """
        ...

    @abstractmethod
    def publish(self, records: List[Record]) -> None:
        """即发即忘：持久化 + 推上传队列，不等待确认。用于高频数据。"""
        ...

    @abstractmethod
    def fork(self) -> "CoreProtocol":
        """
        创建一个新的 CoreProtocol，此函数用于处理多进程fork的情况
        """
        ...

    @abstractmethod
    def deliver_run_finish(self, finish_record: FinishRecord) -> FinishResponse:
        """
        交付运行结束事件，同步事件，等待确认
        """
        ...
