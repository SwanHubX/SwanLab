"""
@author: cunyue
@file: memory.py
@time: 2026/3/31 01:50
@description: 内存信息模块
"""

import subprocess
import sys
from typing import Optional, Tuple

import psutil

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.system import MemorySnapshot, SystemScalar, SystemScalars, SystemShim
from swanlab.sdk.typings.run.system.hardware_vendor import MemoryProtocol
from swanlab.sdk.utils.helper import catch_and_return_none
from swanlab.utils import generate_color


def _bytes_to_snapshot(total_bytes: int, MB: int = 1024**2, GB: int = 1024**3) -> Optional[MemorySnapshot]:
    """将字节数转换为合适单位的 MemorySnapshot，优先 GB，不足 1 GB 时用 MB"""
    gb = total_bytes // GB
    if gb > 0:
        return MemorySnapshot(total=gb, total_unit="GB")
    mb = total_bytes // MB
    if mb > 0:
        return MemorySnapshot(total=mb, total_unit="MB")
    return None


class Memory(MemoryProtocol):
    """
    内存信息模块
    这是通用的内存信息采集实现，全部使用Python标准库，设计上这属于兜底的内存采集实现
    """

    @classmethod
    def new(cls, shim: SystemShim) -> Optional[Tuple["Memory", SystemScalars]]:
        if not shim.enable_memory:
            return None
        if shim.slug == "macos-arm":
            return None
        self = cls(shim)
        scalars: SystemScalars = []
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
        current_process = psutil.Process()
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
    @catch_and_return_none(on_error=lambda e: console.debug("Failed to get memory info: {}", e))
    def get() -> Optional[MemorySnapshot]:
        total_bytes = Memory._get_total_bytes()
        if total_bytes is None or total_bytes <= 0:
            return None
        return _bytes_to_snapshot(total_bytes)

    @staticmethod
    @catch_and_return_none(on_error=lambda e: console.debug("Failed to get total memory bytes: {}", e))
    def _get_total_bytes(timeout_seconds=2.0) -> Optional[int]:
        """获取系统总内存字节数"""
        if sys.platform == "linux":
            with open("/proc/meminfo", "r") as f:
                for line in f:
                    if line.startswith("MemTotal:"):
                        # 格式：MemTotal:       16384000 kB
                        kb = int(line.split()[1])
                        return kb * 1024
        elif sys.platform == "darwin":
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, timeout=timeout_seconds
            )
            if result.returncode == 0:
                return int(result.stdout.strip())
        elif sys.platform == "win32":
            result = subprocess.run(
                ["wmic", "ComputerSystem", "get", "TotalPhysicalMemory"],
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
            )
            if result.returncode == 0:
                lines = [line.strip() for line in result.stdout.strip().split("\n")[1:] if line.strip()]
                if lines:
                    return int(lines[0])

        return None
