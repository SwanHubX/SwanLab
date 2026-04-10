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
from typing import Optional, Tuple

import psutil

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.system import CPUSnapshot, SystemScalar, SystemScalars, SystemShim
from swanlab.sdk.typings.run.system.hardware_vendor import CpuProtocol
from swanlab.sdk.utils.helper import catch_and_return_none
from swanlab.utils import generate_color


class CPU(CpuProtocol):
    """
    CPU 信息模块
    这是通用的 CPU 信息采集实现，全部使用python标准库，设计上这属于兜底的CPU采集实现
    """

    @classmethod
    def new(cls, shim: SystemShim) -> Optional[Tuple["CPU", SystemScalars]]:
        if not shim.enable_cpu:
            return None
        if shim.slug == "macos-arm":
            return None
        self = cls(shim)
        psutil.cpu_percent(interval=None)  # 预热，使后续非阻塞调用能返回真实值
        # 定义所有会被采集的指标
        # 当前进程线程数
        scalars: SystemScalars = []
        # 平均使用率
        usage = SystemScalar(
            key="cpu.pct",
            name="CPU Utilization (%)",
            chart_name="CPU Utilization",
            color=generate_color(0),
        )
        scalars.append(usage)
        self._handlers.append(("cpu.pct", lambda: psutil.cpu_percent(interval=None)))
        # 当前线程数
        threads = SystemScalar(
            key="cpu.thds",
            name="Process CPU Threads",
            chart_name="Process CPU Threads",
            color=generate_color(0),
        )
        scalars.append(threads)
        self._handlers.append(("cpu.thds", lambda: psutil.Process().num_threads()))
        return self, scalars

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
