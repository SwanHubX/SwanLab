"""
@author: cunyue
@file: apple.py
@time: 2026/3/31 01:51
@description: 苹果统一芯片（Apple Silicon）相关的硬件供应商信息和功能实现
对于原本的Intel架构的Mac，我们在cpu.py和memory.py中统一采集信息
"""

import json
import multiprocessing
import subprocess
import sys
from typing import Optional, Tuple

import psutil

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.system import AppleSiliconSnapshot, SystemScalar, SystemScalars, SystemShim
from swanlab.sdk.typings.run.system.hardware_vendor import AppleSiliconProtocol
from swanlab.sdk.utils.helper import catch_and_return_none
from swanlab.utils import generate_color


class Apple(AppleSiliconProtocol):
    """
    苹果统一芯片（Apple Silicon）相关的硬件供应商信息和功能实现
    在macOS-arm平台上替代CPU和Memory模块，统一采集CPU、内存相关指标
    """

    @classmethod
    def new(cls, shim: SystemShim) -> Optional[Tuple["Apple", SystemScalars]]:
        # 只有macOS-arm架构才支持苹果统一芯片
        if shim.slug != "macos-arm":
            return None
        self = cls(shim)
        psutil.cpu_percent(interval=None)  # 预热，使后续非阻塞调用能返回真实值
        scalars: SystemScalars = []
        current_process = psutil.Process()
        if shim.enable_cpu:
            # CPU 使用率
            cpu_pct = SystemScalar(
                key="cpu.pct",
                name="CPU Utilization (%)",
                chart_name="CPU Utilization",
                color=generate_color(0),
            )
            scalars.append(cpu_pct)
            self._handlers.append(("cpu.pct", lambda: psutil.cpu_percent(interval=None)))
            # 当前进程线程数
            cpu_thds = SystemScalar(
                key="cpu.thds",
                name="Process CPU Threads",
                chart_name="Process CPU Threads",
                color=generate_color(0),
            )
            scalars.append(cpu_thds)
            self._handlers.append(("cpu.thds", lambda: current_process.num_threads()))
        if shim.enable_memory:
            # 系统内存使用率
            mem_pct = SystemScalar(
                key="mem.pct",
                name="System Memory Utilization (%)",
                chart_name="System Memory",
                color=generate_color(0),
            )
            scalars.append(mem_pct)
            self._handlers.append(("mem.pct", lambda: psutil.virtual_memory().percent))
            # 当前进程内存使用（RSS，非交换区）
            mem_proc = SystemScalar(
                key="mem.proc",
                name="Process Memory In Use (non-swap) (MB)",
                chart_name="Process Memory",
                color=generate_color(0),
            )
            scalars.append(mem_proc)
            self._handlers.append(("mem.proc", lambda: current_process.memory_info().rss / 1024 / 1024))
        return self, scalars

    @staticmethod
    @catch_and_return_none(on_error=lambda e: console.debug("Failed to get Apple Silicon info: {}", e))
    def get() -> Optional[AppleSiliconSnapshot]:
        if sys.platform != "darwin":
            return None
        result = subprocess.run(
            ["system_profiler", "SPHardwareDataType", "-json"], capture_output=True, text=True, timeout=5
        )
        hardware_info = json.loads(result.stdout)["SPHardwareDataType"][0]
        chip_type = hardware_info.get("chip_type") or hardware_info.get("cpu_type")
        if chip_type is None:
            return None
        memory = str(hardware_info["physical_memory"]).lower().replace("gb", "").strip()
        cpu_count = multiprocessing.cpu_count()
        return AppleSiliconSnapshot(name=chip_type, memory=int(memory), memory_unit="GB", cpu_count=cpu_count)
