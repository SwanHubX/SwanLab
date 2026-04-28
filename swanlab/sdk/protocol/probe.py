"""
@author: cunyue
@file: probe.py
@time: 2026/4/21 15:13
@description: SwanLab Probe 接口协议，负责启动系统信息采集、监控服务
"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import TYPE_CHECKING

from swanlab.sdk.internal.pkg import safe

if TYPE_CHECKING:
    from swanlab.sdk.internal.context import RunContext


class ProbeEnum(str, Enum):
    PROBE_PYTHON = "ProbePython"
    PROBE = "Probe"


class ProbeProtocol(ABC):
    """
    SwanLab Probe 协议
    协议层作为一个“垫片”，抹平前端与后端的实现爱你差异
    Probe设计上作为一个微服务，仅作为后台进程启动，本身与core直接交互而非SDK本体，本体的作用是启动、关闭probe
    """

    def __init__(self, ctx: "RunContext"):
        self._ctx = ctx
        self._mode = ctx.config.settings.mode

    def start(self):
        """
        启动 probe 服务，首先做相关上下文检查，然后启动服务
        """
        # 检查 probe 所需的上下文信息
        assert self._ctx.core is not None, "core must be started before starting probe"
        if self._mode == "disabled":
            return self._start_when_disabled()
        assert self._ctx.files_dir.is_dir(), "files_dir must be set before starting probe"
        with safe.block(message="probe start error"):
            if self._mode == "online":
                return self._start_when_online()
            elif self._mode == "local":
                return self._start_when_local()
            return self._start_when_offline()

    def _start_when_disabled(self) -> None: ...

    @abstractmethod
    def _start_when_local(self) -> None: ...

    @abstractmethod
    def _start_when_offline(self) -> None: ...

    @abstractmethod
    def _start_when_online(self) -> None: ...

    def finish(self):
        """
        停止 probe 服务
        """
        with safe.block(message="probe finish error"):
            if self._mode == "online":
                return self._finish_when_online()
            elif self._mode == "local":
                return self._finish_when_local()
            elif self._mode == "offline":
                return self._finish_when_offline()
            return self._finish_when_disabled()

    def _finish_when_disabled(self) -> None: ...

    @abstractmethod
    def _finish_when_local(self) -> None: ...

    @abstractmethod
    def _finish_when_offline(self) -> None: ...

    @abstractmethod
    def _finish_when_online(self) -> None: ...
