"""
@author: cunyue
@file: probe.py
@time: 2026/4/21 15:13
@description: SwanLab Probe 接口协议，负责启动系统信息采集、监控服务
"""

from abc import ABC, abstractmethod
from enum import Enum

from swanlab.proto.swanlab.grpc.probe.v1.probe_pb2 import DeliverProbeStartRequest, GetMetadataSnapshotResponse
from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.typings.run import ModeType


class ProbeEnum(str, Enum):
    PROBE_PYTHON = "ProbePython"
    PROBE = "Probe"


class ProbeProtocol(ABC):
    """
    SwanLab Probe 协议
    协议层作为一个“垫片”，抹平前端与后端的实现差异
    Probe设计上作为一个微服务，仅作为后台进程启动，本身与core直接交互而非SDK本体，本体的作用是启动、关闭probe
    """

    def __init__(self, mode: ModeType):
        self._mode = mode

    def deliver_probe_start(self, start_request: DeliverProbeStartRequest) -> None:
        """
        启动 probe 服务，首先做相关上下文检查，然后启动服务
        """
        if self._mode == "disabled":
            return self._start_when_disabled()
        with safe.block(message="probe start error"):
            return self._start_when_enabled(start_request)

    def _start_when_disabled(self) -> None: ...

    @abstractmethod
    def _start_when_enabled(self, start_request: DeliverProbeStartRequest) -> None: ...

    def get_metadata_snapshot(self) -> GetMetadataSnapshotResponse:
        if self._mode == "disabled":
            return self._get_metadata_snapshot_when_disabled()
        with safe.block(message="get metadata snapshot error"):
            return self._get_metadata_snapshot_when_enabled()
        return GetMetadataSnapshotResponse(success=False, message="Unknown error.")

    def _get_metadata_snapshot_when_disabled(self) -> GetMetadataSnapshotResponse:
        return GetMetadataSnapshotResponse(success=False, message="I'm a teapot.")

    @abstractmethod
    def _get_metadata_snapshot_when_enabled(self) -> GetMetadataSnapshotResponse: ...

    def deliver_probe_finish(self) -> None:
        """
        停止 probe 服务
        """
        if self._mode == "disabled":
            return self._finish_when_disabled()
        with safe.block(message="probe finish error"):
            return self._finish_when_enabled()

    def _finish_when_disabled(self) -> None: ...

    @abstractmethod
    def _finish_when_enabled(self) -> None: ...
