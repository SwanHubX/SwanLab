"""
@author: cunyue
@file: utils.py
@time: 2024/12/9 16:33
@description: 硬件信息采集工具函数
"""

import random

import psutil
from swankit.callback.models import ColumnConfig

from .type import HardwareInfo, HardwareInfoList

ALPHABET = "abcdefghijklmnopqrstuvwxyz0123456789"


def random_index(length: int = 8) -> str:
    """
    随机生成八位字符串，用于标识图表的index
    """
    return "".join(random.choices(ALPHABET, k=length))


def generate_key(suffix: str, length: int = 4) -> str:
    """
    生成key，用于标识系统列，避免与用户输入的key冲突
    """
    return "".join(random.choices(ALPHABET, k=length)) + "." + suffix


class CpuCollector:
    """
    cpu采集基类，为子类赋予cpu采集的能力
    """

    CPU_CONFIG = ColumnConfig(y_range=(0, 100), chart_name="CPU Utilization (%)")
    PER_CPU_CONFIG = ColumnConfig(y_range=(0, 100), chart_name="CPU Utilization (per core) (%)")
    THDS_CONFIG = ColumnConfig(y_range=(0, None), chart_name="Process CPU Threads")

    per_cpu_configs = []
    # 随机生成一个index，用于标识图表
    per_cpu_usage_chart_index = random_index()
    cpu_usage_key = generate_key("cpu.pct")
    per_cpu_usage_key = generate_key("cpu.{idx}.pct")
    proc_thds_key = generate_key("cpu.thds")

    def get_cpu_usage(self) -> HardwareInfo:
        """
        获取当前 CPU 使用率
        """
        return {
            "key": self.cpu_usage_key,
            "name": "CPU Utilization (%)",
            "value": psutil.cpu_percent(interval=1),
            "config": self.CPU_CONFIG,
        }

    def get_per_cpu_usage(self) -> HardwareInfoList:
        """
        获取每个 CPU 核心的使用率
        """
        per_cpu_usage = psutil.cpu_percent(interval=1, percpu=True)
        result: HardwareInfoList = []
        # 避免每次调用都创建新的配置
        if len(self.per_cpu_configs) != len(per_cpu_usage):
            self.per_cpu_configs = [
                self.PER_CPU_CONFIG.clone(metric_name=f"CPU {idx}", chart_index=self.per_cpu_usage_chart_index)
                for idx in range(len(per_cpu_usage))
            ]
        for idx, value in enumerate(per_cpu_usage):
            info: HardwareInfo = {
                "key": self.per_cpu_usage_key.format(idx=idx),
                "name": f"CPU {idx} Utilization (%)",
                "value": value,
                "config": self.per_cpu_configs[idx],
            }
            result.append(info)
        return result

    def get_cur_proc_thds_num(self, proc: psutil.Process) -> HardwareInfo:
        """
        获取当前进程的线程数
        """
        return {
            "key": self.proc_thds_key,
            "name": "Process CPU Threads",
            "value": proc.num_threads(),
            "config": self.THDS_CONFIG,
        }


class MemoryCollector:
    """
    内存采集基类，为子类赋予内存采集的能力
    """

    MB = 1024 * 1024
    MEM_CONFIG = ColumnConfig(y_range=(0, 100))
    PROC_MEM_PCT_CONFIG = ColumnConfig(y_range=(0, 100))

    mem_usage_key = generate_key("mem.pct")
    mem_proc_key = generate_key("mem.proc")
    mem_proc_pct_key = generate_key("mem.proc.pct")
    mem_proc_avail_key = generate_key("mem.proc.avail")

    def get_mem_usage(self) -> HardwareInfo:
        """
        获取当前系统内存使用率
        """
        return {
            "key": self.mem_usage_key,
            "name": "System Memory Utilization (%)",
            "value": psutil.virtual_memory().percent,
            "config": self.MEM_CONFIG,
        }

    def get_cur_proc_mem(self, proc: psutil.Process) -> HardwareInfoList:
        """
        获取当前进程的内存使用情况
        """
        mem_info = proc.memory_info()
        virtual_memory = psutil.virtual_memory()
        mem_proc: HardwareInfo = {
            "key": self.mem_proc_key,
            "name": "Process Memory In Use (non-swap) (MB)",
            "value": mem_info.rss / self.MB,
            "config": None,
        }
        mem_proc_pct: HardwareInfo = {
            "key": self.mem_proc_pct_key,
            "name": "Process Memory Utilization (%)",
            "value": proc.memory_percent(),
            "config": self.PROC_MEM_PCT_CONFIG,
        }
        mem_proc_avail: HardwareInfo = {
            "key": self.mem_proc_avail_key,
            "name": "Process Memory Available (MB)",
            "value": virtual_memory.available / self.MB,
            "config": None,
        }

        return [mem_proc, mem_proc_pct, mem_proc_avail]
