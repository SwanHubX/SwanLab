"""
@author: cunyue
@file: cpu.py
@time: 2026/3/31 01:50
@description: CPU 信息模块
"""

import multiprocessing
import platform
import subprocess
import sys
from typing import TYPE_CHECKING, Optional

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.system import CPUSnapshot, SystemShim
from swanlab.sdk.typings.run.system.hardware_vendor import CpuProtocol
from swanlab.sdk.utils.helper import catch_and_return_none

if TYPE_CHECKING:
    from swanlab import Run


class CPU(CpuProtocol):
    """
    CPU 信息模块
    这是通用的 CPU 信息采集实现，全部使用python标准库，设计上这属于兜底的CPU采集实现
    """

    @classmethod
    def new(cls, run: "Run", shim: SystemShim) -> Optional["CPU"]:
        if shim.slug == "macos-arm":
            return None
        return cls(shim)

    def collect(self):
        raise NotImplementedError()

    @staticmethod
    @catch_and_return_none(on_error=lambda e: console.debug("Failed to get CPU info: {}", e))
    def get() -> Optional[CPUSnapshot]:
        # 1. 获取 CPU 品牌
        brand = CPU._get_real_brand() or platform.processor() or None
        # 2. 获取 CPU 核心数
        logical_count = multiprocessing.cpu_count()

        # 尝试获取物理核心数
        physical_count = None
        timeout_seconds = 2.0  # 防止子进程卡死主线程

        if sys.platform == "linux":
            # 引入 Socket 概念防止多路 CPU Core ID 重复导致计算偏小
            result = subprocess.run(
                ["lscpu", "-p=Core,Socket"], capture_output=True, text=True, timeout=timeout_seconds
            )
            if result.returncode == 0:
                physical_count = len(
                    set(line for line in result.stdout.strip().split("\n") if not line.startswith("#") and line.strip())
                )
        elif sys.platform == "darwin":
            result = subprocess.run(
                ["sysctl", "-n", "hw.physicalcpu"], capture_output=True, text=True, timeout=timeout_seconds
            )
            if result.returncode == 0:
                physical_count = int(result.stdout.strip())
        elif sys.platform == "win32":
            result = subprocess.run(
                ["wmic", "cpu", "get", "NumberOfCores"], capture_output=True, text=True, timeout=timeout_seconds
            )
            if result.returncode == 0:
                cores = [line.strip() for line in result.stdout.strip().split("\n")[1:] if line.strip()]
                physical_count = sum(int(c) for c in cores if c.isdigit())

        return CPUSnapshot(brand=brand, physical_count=physical_count, logical_count=logical_count)

    @staticmethod
    @catch_and_return_none(on_error=lambda e: console.debug("Failed to get real CPU brand: {}", e))
    def _get_real_brand() -> Optional[str]:
        """尝试获取真实的 CPU 品牌名称"""
        if sys.platform == "linux":
            with open("/proc/cpuinfo", "r") as f:
                for line in f:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
        elif sys.platform == "darwin":
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, timeout=2.0
            )
            if result.returncode == 0:
                return result.stdout.strip()
        elif sys.platform == "win32":
            import winreg

            key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
            brand, _ = winreg.QueryValueEx(key, "ProcessorNameString")
            return brand.strip()
