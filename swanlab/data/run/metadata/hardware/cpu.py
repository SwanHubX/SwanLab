"""
@author: cunyue
@file: cpu.py
@time: 2024/12/3 20:12
@description: CPU信息采集
"""

import multiprocessing
import platform
import subprocess

import psutil
import os
from swankit.env import is_macos

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult, HardwareCollector, HardwareInfoList
from .utils import CpuCollector as C


def get_cpu_info() -> HardwareFuncResult:
    """获取 CPU 信息"""
    if is_macos():
        return None, None
    info = {"brand": None, "cores": None}

    # 获取 CPU 品牌, 根据不同操作系统调用不同的函数
    if platform.system() == "Windows":
        info["brand"] = get_cpu_brand_windows()
    elif platform.system() == "Linux":
        info["brand"] = get_cpu_brand_linux()
    else:
        # 其他情况，暂时不支持
        return None, None
    try:
        # 获取 CPU 核心数
        info["cores"] = multiprocessing.cpu_count()
    except Exception:  # noqa
        pass

    return info, CpuCollector()


def get_cpu_brand_windows():
    try:
        # 使用 WMIC 命令获取 CPU 品牌
        result = subprocess.run(["wmic", "cpu", "get", "name"], capture_output=True, text=True)
        cpu_brand = result.stdout.strip().split("\n")[-1].strip()
        return cpu_brand
    except Exception:  # noqa
        return None


def get_cpu_brand_linux():
    # lscpu 命令获取 CPU 品牌
    cpu_brand = None
    try:
        result = subprocess.run(["lscpu"], capture_output=True, text=True)
        for line in result.stdout.split("\n"):
            if "model name" in line.lower():
                cpu_brand = line.split(":")[1].strip()
                break
    except Exception:  # noqa
        pass
    if cpu_brand is None:
        try:
            # 使用 cat /proc/cpuinfo 获取CPU品牌
            with open("/proc/cpuinfo", "r") as f:
                for line in f:
                    if "model name" in line.lower():
                        cpu_brand = line.split(":")[1].strip()
                        break
        except Exception:  # noqa
            pass
    return cpu_brand


class CpuCollector(HardwareCollector, C):
    def __init__(self):
        super().__init__()
        self.current_process = psutil.Process()

    def collect(self) -> HardwareInfoList:
        return [
            self.get_cpu_usage(),
            # *self.get_per_cpu_usage(),
            self.get_cur_proc_thds_num(self.current_process),
        ]
