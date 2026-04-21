"""
@author: cunyue
@file: __init__.py.py
@time: 2026/4/21 15:11
@description: SwanLab SDK 内部系统组件，如硬件监控、元数据采集、终端日志收集等
与core设计理念一致，未来作为单独微服务实现
"""

from typing import Optional

from swanlab.sdk.internal.context import RunContext
from swanlab.sdk.internal.probe_python.monitor import Monitor
from swanlab.sdk.protocol.probe import ProbeProtocol


class ProbePython(ProbeProtocol):
    def __init__(self, ctx: RunContext):
        super().__init__(ctx)
        self._monitor: Optional[Monitor] = None

    def _start(self):
        """
        启动硬件采集
        """
        # settings = self._ctx.config.settings
        # # 1. 系统环境信息采集
        # git_snapshot = git.get() if settings.environment.git else None
        # runtime_snapshot = runtime.get() if settings.environment.runtime else None
        # conda_snapshot = conda.get() if settings.environment.conda else None
        # requirements_snapshot = requirements.get() if settings.environment.requirements else None
        # # 2. 系统硬件信息采集，只有在硬件采集和监控都不开启时才不采集
        # hardware_snapshot = None
        # if settings.environment.hardware or settings.monitor.enable:
        #     # 2.1 苹果
        #     apple_snapshot = Apple.get()
        #     # 2.2 通用硬件
        #     cpu_snapshot = None
        #     memory_snapshot = None
        #     if apple_snapshot is None:
        #         cpu_snapshot = CPU.get()
        #         memory_snapshot = Memory.get()
        #     # 2.3 各家加速器厂商
        #
        #     # 3. 组合数据
        #     hardware_snapshot = HardwareSnapshot(
        #         cpu=cpu_snapshot,
        #         memory=memory_snapshot,
        #         apple_silicon=apple_snapshot,
        #     )
        # metadata = MetadataSnapshot(
        #     hardware=hardware_snapshot,
        #     git=git_snapshot,
        #     runtime=runtime_snapshot,
        # )
        # system_shim = SystemShim.from_snapshot(metadata, platform=sys.platform)
        # # 如果硬采集被设置为False，删除硬件信息
        # if not settings.environment.hardware:
        #     metadata = metadata.del_hardware()
        #
        # system_environment = SystemEnvironment(
        #     shim=system_shim,
        #     metadata=metadata,
        #     requirements=requirements_snapshot,
        #     conda=conda_snapshot,
        # )
        # # 3. 硬件监控
        # if settings.monitor.enable:
        #     self._monitor = Monitor(system_shim)

    def _start_when_local(self) -> None:
        pass

    def _start_when_offline(self) -> None:
        pass

    def _start_when_cloud(self) -> None:
        pass

    def _finish_when_local(self) -> None:
        pass

    def _finish_when_offline(self) -> None:
        pass

    def _finish_when_cloud(self) -> None:
        pass


# def create_monitor(ctx: RunContext, e: EmitterProtocol):
#     if ctx.config.settings.mode == "disabled":
#         return None
#     this_monitor: Optional["Monitor"] = None
#     sys_info, hardware_monitor = _new_raw(ctx)
#     ts = Timestamp()
#     ts.GetCurrentTime()
#     if sys_info.metadata:
#         fs.safe_write(ctx.metadata_file, sys_info.metadata.model_dump_json())
#         e.emit(MetadataEvent(timestamp=ts))
#     if sys_info.requirements:
#         fs.safe_write(ctx.requirements_file, sys_info.requirements)
#         e.emit(RequirementsEvent(timestamp=ts))
#     if sys_info.conda:
#         fs.safe_write(ctx.conda_file, sys_info.conda)
#         e.emit(CondaEvent(timestamp=ts))
#     if hardware_monitor is not None and hardware_monitor.start(ctx, e):
#         this_monitor = hardware_monitor
#     return this_monitor
