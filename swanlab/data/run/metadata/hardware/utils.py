"""
@author: cunyue
@file: utils.py
@time: 2024/12/9 16:33
@description: 硬件信息采集工具函数
"""

import random

import psutil

from .type import HardwareInfo, HardwareInfoList, HardwareConfig

ALPHABET = "abcdefghijklmnopqrstuvwxyz0123456789"


def random_index(length: int = 8) -> str:
    """
    随机生成八位字符串，用于标识图表的index
    """
    return "".join(random.choices(ALPHABET, k=length))


def generate_key(suffix: str) -> str:
    """
    生成key，用于标识系统列，避免与用户输入的key冲突
    """
    return "__swanlab__." + suffix


# CPU 使用率
CPU_PCT_KEY = generate_key("cpu.pct")
CPU_PCT_CONFIG = HardwareConfig(
    y_range=(0, 100),
    chart_name="CPU Utilization (%)",
).clone()

# CPU 核心使用率
CPU_INDEX_PER_KEY = generate_key("cpu.{idx}.pct")
CPU_INDEX_PER_CONFIG = HardwareConfig(
    y_range=(0, 100),
    chart_name="CPU Utilization (per core) (%)",
    chart_index=random_index(),
)
CPU_INDEX_PER_CONFIGS = [CPU_INDEX_PER_CONFIG.clone(metric_name=f"CPU {idx}") for idx in range(psutil.cpu_count())]

# CPU 线程数
CPU_THDS_KEY = generate_key("cpu.thds")
CPU_THDS_CONFIG = HardwareConfig(
    y_range=(0, None),
    chart_name="Process CPU Threads",
).clone()


class CpuCollector:
    """
    cpu采集基类，为子类赋予cpu采集的能力
    """

    __per_cpu_configs = []

    @staticmethod
    def get_cpu_usage() -> HardwareInfo:
        """
        获取当前 CPU 使用率
        """
        return {
            "key": CPU_PCT_KEY,
            "name": "CPU Utilization (%)",
            "value": psutil.cpu_percent(interval=1),
            "config": CPU_PCT_CONFIG,
        }

    @staticmethod
    def get_per_cpu_usage() -> HardwareInfoList:
        """
        获取每个 CPU 核心的使用率
        """
        per_cpu_usages = psutil.cpu_percent(interval=1, percpu=True)
        result: HardwareInfoList = []
        for idx, value in enumerate(per_cpu_usages):
            info: HardwareInfo = {
                "key": CPU_INDEX_PER_KEY.format(idx=idx),
                "name": f"CPU {idx} Utilization (%)",
                "value": value,
                "config": CPU_INDEX_PER_CONFIGS[idx],
            }
            result.append(info)
        return result

    @staticmethod
    def get_cur_proc_thds_num(proc: psutil.Process) -> HardwareInfo:
        """
        获取当前进程的线程数
        """
        return {
            "key": CPU_THDS_KEY,
            "name": "Process CPU Threads",
            "value": proc.num_threads(),
            "config": CPU_THDS_CONFIG,
        }


# 内存使用率
MEM_PCT_KEY = generate_key("mem.pct")
MEM_PCT_CONFIG = HardwareConfig(
    y_range=(0, 100),
).clone()

# 进程内存使用情况
MEM_PROC_KEY = generate_key("mem.proc")
MEM_PROC_CONFIG = HardwareConfig(
    y_range=(0, None),
    chart_name="Process Memory In Use (non-swap) (MB)",
).clone()

# 进程内存使用率
MEM_PROC_PCT_KEY = generate_key("mem.proc.pct")
MEM_PROC_PCT_CONFIG = HardwareConfig(
    y_range=(0, 100),
    chart_name="Process Memory Utilization (%)",
).clone()

# 进程内存可用情况
PROC_MEM_AVAIL_KEY = generate_key("mem.proc.avail")
PROC_MEM_AVAIL_CONFIG = HardwareConfig(
    y_range=(0, None),
    chart_name="Process Memory Available (non-swap) (MB)",
).clone()


class MemoryCollector:
    """
    内存采集基类，为子类赋予内存采集的能力
    """

    @staticmethod
    def get_mem_usage() -> HardwareInfo:
        """
        获取当前系统内存使用率
        """
        return {
            "key": MEM_PCT_KEY,
            "name": "System Memory Utilization (%)",
            "value": psutil.virtual_memory().percent,
            "config": MEM_PCT_CONFIG,
        }

    @staticmethod
    def get_cur_proc_mem(proc: psutil.Process) -> HardwareInfoList:
        """
        获取当前进程的内存使用情况
        """
        mem_info = proc.memory_info()
        virtual_memory = psutil.virtual_memory()
        mem_proc: HardwareInfo = {
            "key": MEM_PROC_KEY,
            "name": "Process Memory In Use (non-swap) (MB)",
            "value": mem_info.rss / 1024 / 1024,
            "config": MEM_PROC_CONFIG,
        }
        mem_proc_pct: HardwareInfo = {
            "key": MEM_PROC_PCT_KEY,
            "name": "Process Memory Utilization (%)",
            "value": proc.memory_percent(),
            "config": MEM_PROC_PCT_CONFIG,
        }
        mem_proc_avail: HardwareInfo = {
            "key": PROC_MEM_AVAIL_KEY,
            "name": "Process Memory Available (non-swap) (MB)",
            "value": virtual_memory.available / 1024 / 1024,
            "config": PROC_MEM_AVAIL_CONFIG,
        }

        return [mem_proc, mem_proc_pct, mem_proc_avail]
