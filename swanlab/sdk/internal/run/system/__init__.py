"""
@author: cunyue
@file: __init__.py.py
@time: 2026/3/11 17:55
@description: SwanLab SDK 内部系统组件，如硬件监控、元数据采集、终端日志收集等
"""

import sys
from typing import Optional

from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.sdk.internal.bus import CondaEvent, MetadataEvent, RequirementsEvent
from swanlab.sdk.internal.bus.emitter import EmitterProtocol
from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.internal.run.system.environment import conda, git, requirements, runtime
from swanlab.sdk.internal.run.system.hardware_vendor.apple import Apple
from swanlab.sdk.internal.run.system.hardware_vendor.cpu import CPU
from swanlab.sdk.internal.run.system.hardware_vendor.memory import Memory
from swanlab.sdk.internal.run.system.monitor import Monitor
from swanlab.sdk.typings.run.system import HardwareSnapshot, MetadataSnapshot, SystemEnvironment, SystemShim

__all__ = ["create_monitor", "Monitor"]


def _new_raw(ctx: RunContext):
    """
    创建硬件监控模块
    """
    settings = ctx.config.settings
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

        # 3. 组合数据
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

    system_environment = SystemEnvironment(
        shim=system_shim,
        metadata=metadata,
        requirements=requirements_snapshot,
        conda=conda_snapshot,
    )
    # 3. 硬件监控
    hardware_monitor: Optional["Monitor"] = None
    if settings.monitor.enable:
        hardware_monitor = Monitor(system_shim)
    return system_environment, hardware_monitor


def create_monitor(ctx: RunContext, e: EmitterProtocol):
    if ctx.config.settings.mode == "disabled":
        return None
    this_monitor: Optional["Monitor"] = None
    sys_info, hardware_monitor = _new_raw(ctx)
    ts = Timestamp()
    ts.GetCurrentTime()
    if sys_info.metadata:
        fs.safe_write(ctx.metadata_file, sys_info.metadata.model_dump_json())
        e.emit(MetadataEvent(timestamp=ts))
    if sys_info.requirements:
        fs.safe_write(ctx.requirements_file, sys_info.requirements)
        e.emit(RequirementsEvent(timestamp=ts))
    if sys_info.conda:
        fs.safe_write(ctx.conda_file, sys_info.conda)
        e.emit(CondaEvent(timestamp=ts))
    if hardware_monitor is not None and hardware_monitor.start(ctx, e):
        this_monitor = hardware_monitor
    return this_monitor
