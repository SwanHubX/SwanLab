"""
@author: cunyue
@file: __init__.py.py
@time: 2026/4/21 15:11
@description: SwanLab SDK 内部系统组件，如硬件监控、元数据采集、终端日志收集等
与core设计理念一致，未来作为单独微服务实现
"""

import sys
from typing import Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.system.v1.env_pb2 import CondaRecord, MetadataRecord, RequirementsRecord
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.internal.probe_python.environment import conda, git, requirements, runtime
from swanlab.sdk.internal.probe_python.hardware_vendor.apple import Apple
from swanlab.sdk.internal.probe_python.hardware_vendor.cpu import CPU
from swanlab.sdk.internal.probe_python.hardware_vendor.memory import Memory
from swanlab.sdk.internal.probe_python.monitor import Monitor
from swanlab.sdk.protocol.probe import ProbeProtocol
from swanlab.sdk.typings.probe_python import HardwareSnapshot, MetadataSnapshot, SystemEnvironment, SystemShim


class ProbePython(ProbeProtocol):
    def __init__(self, ctx: RunContext):
        super().__init__(ctx)
        # self._monitor is not None 作为硬件监控是否开启的标记
        self._monitor: Optional[Monitor] = None

    def _start(self):
        """
        启动硬件采集
        """
        settings = self._ctx.config.settings
        # 1. 系统环境信息采集
        git_snapshot = git.get() if settings.environment.git else None
        runtime_snapshot = runtime.get() if settings.environment.runtime else None
        conda_snapshot = conda.get() if settings.environment.conda else None
        requirements_snapshot = requirements.get() if settings.environment.requirements else None
        # 2. 系统硬件信息采集，只有在硬件采集和监控都不开启时才不采集
        hardware_snapshot = None
        if settings.environment.hardware or settings.monitor.enable:
            # 2.1 苹果
            apple_snapshot = Apple.get()
            # 2.2 通用硬件
            cpu_snapshot = None
            memory_snapshot = None
            if apple_snapshot is None:
                cpu_snapshot = CPU.get()
                memory_snapshot = Memory.get()
            # 2.3 各家加速器厂商

            # 2.4 组合数据
            hardware_snapshot = HardwareSnapshot(
                cpu=cpu_snapshot,
                memory=memory_snapshot,
                apple_silicon=apple_snapshot,
            )
        metadata = MetadataSnapshot(
            hardware=hardware_snapshot,
            git=git_snapshot,
            runtime=runtime_snapshot,
        )
        system_shim = SystemShim.from_snapshot(metadata, platform=sys.platform)
        # 如果硬采集被设置为False，删除硬件信息
        if not settings.environment.hardware:
            metadata = metadata.del_hardware()

        sys_info = SystemEnvironment(
            shim=system_shim,
            metadata=metadata,
            requirements=requirements_snapshot,
            conda=conda_snapshot,
        )
        # 4. 向core发送记录
        ts = Timestamp()
        ts.GetCurrentTime()
        if sys_info.metadata:
            fs.safe_write(self._ctx.metadata_file, sys_info.metadata.model_dump_json(by_alias=True))
            self._ctx.core.upsert_metadata([MetadataRecord(timestamp=ts)])
        if sys_info.requirements:
            fs.safe_write(self._ctx.requirements_file, sys_info.requirements)
            self._ctx.core.upsert_requirements([RequirementsRecord(timestamp=ts)])
        if sys_info.conda:
            fs.safe_write(self._ctx.conda_file, sys_info.conda)
            self._ctx.core.upsert_conda([CondaRecord(timestamp=ts)])

        # 5. 设置硬件监控
        if settings.monitor.enable:
            if (hm := Monitor(system_shim, self._ctx.core).start(self._ctx)) is not None:
                self._monitor = hm

    def _start_when_local(self) -> None:
        self._start()

    def _start_when_offline(self) -> None:
        self._start()

    def _start_when_cloud(self) -> None:
        self._start()

    def _finish(self):
        if self._monitor is not None:
            self._monitor.stop()

    def _finish_when_local(self) -> None:
        self._finish()

    def _finish_when_offline(self) -> None:
        if self._monitor is not None:
            self._monitor.stop()

    def _finish_when_cloud(self) -> None:
        if self._monitor is not None:
            self._monitor.stop()
