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
from swanlab.sdk.internal.pkg import safe

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext

__all__ = ["CoreEnum", "CoreProtocol"]


class CoreEnum(str, Enum):
    CORE_PYTHON = "CorePython"
    CORE = "Core"


class CoreProtocol(ABC):
    """
    SwanLab Core 协议
    协议层作为一个"垫片"，抹平前端与后端实现的差异，我们对CoreProtocol的错误处理理念基本与go一致，不raise Error，仅返回处理结果，上层根据返回结果做相应处理
    """

    def __init__(self, ctx: "RunContext"):
        self._ctx = ctx
        self._mode = ctx.config.settings.mode

    def deliver_run_start(self, start_record: StartRecord) -> StartResponse:
        """
        交付运行开始事件，同步事件，等待确认
        约定应该在此处完成Core相关组件的初始化，在这里完成不同模式的初始化函数分发
        """
        with safe.block(message="run start error"):
            if self._mode == "cloud":
                return self._start_when_cloud(start_record)
            elif self._mode == "local":
                return self._start_when_local(start_record)
            elif self._mode == "offline":
                return self._start_when_offline(start_record)
            return self._start_when_disabled(start_record)
        return StartResponse(success=False, message="Failed to start run")

    @staticmethod
    def _start_when_disabled(start_record: StartRecord) -> StartResponse:
        return StartResponse(success=True, message="I'm a teapot.", run=start_record)

    @abstractmethod
    def _start_when_local(self, start_record: StartRecord) -> StartResponse: ...

    @abstractmethod
    def _start_when_offline(self, start_record: StartRecord) -> StartResponse: ...

    @abstractmethod
    def _start_when_cloud(self, start_record: StartRecord) -> StartResponse: ...

    def publish(self, records: List[Record]) -> None:
        """
        即发即忘：持久化 + 推上传队列，不等待确认。用于高频数据。
        在这里完成不同模式的发布函数分发
        """
        with safe.block(message="publish error"):
            if self._mode == "cloud":
                return self._publish_when_cloud(records)
            elif self._mode == "local":
                return self._publish_when_local(records)
            elif self._mode == "offline":
                return self._publish_when_offline(records)
            return self._publish_when_disabled(records)

    def _publish_when_disabled(self, records: List[Record]) -> None: ...

    @abstractmethod
    def _publish_when_local(self, records: List[Record]) -> None: ...

    @abstractmethod
    def _publish_when_offline(self, records: List[Record]) -> None: ...

    @abstractmethod
    def _publish_when_cloud(self, records: List[Record]) -> None: ...

    @abstractmethod
    def fork(self) -> "CoreProtocol":
        """
        创建一个新的 CoreProtocol，此函数用于处理多进程fork的情况
        """
        ...

    def deliver_run_finish(self, finish_record: FinishRecord) -> FinishResponse:
        """
        交付运行结束事件，同步事件，等待确认
        约定应该在此处完成Core相关组件的清理，在这里完成不同模式的结束函数分发
        """
        with safe.block(message="run finish error"):
            if self._mode == "cloud":
                return self._finish_when_cloud(finish_record)
            elif self._mode == "local":
                return self._finish_when_local(finish_record)
            elif self._mode == "offline":
                return self._finish_when_offline(finish_record)
            return self._finish_when_disabled()
        return FinishResponse(success=False, message="Failed to finish run")

    @staticmethod
    def _finish_when_disabled() -> FinishResponse:
        return FinishResponse(success=True, message="I'm a teapot.")

    @abstractmethod
    def _finish_when_local(self, finish_record: FinishRecord) -> FinishResponse: ...

    @abstractmethod
    def _finish_when_offline(self, finish_record: FinishRecord) -> FinishResponse: ...

    @abstractmethod
    def _finish_when_cloud(self, finish_record: FinishRecord) -> FinishResponse: ...
