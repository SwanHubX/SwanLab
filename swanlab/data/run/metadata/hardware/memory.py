"""
@author: cunyue
@file: memory.py
@time: 2024/12/3 20:13
@description: 内存信息采集
"""

from typing import List

import psutil
from swankit.env import is_macos

from .type import HardwareFuncResult, HardwareCollector, HardwareInfo
from .utils import MemoryCollector as M


def get_memory_size() -> HardwareFuncResult:
    """获取内存大小"""
    if is_macos():
        return None, None
    try:
        # 获取系统总内存大小
        total = psutil.virtual_memory().total
        total_memory = round(total / (1024**3))  # 单位为GB
        return total_memory, MemoryCollector()
    except Exception:  # noqa
        return None, None


class MemoryCollector(HardwareCollector, M):
    def __init__(self):
        super().__init__()
        self.current_process = psutil.Process()

    def collect(self) -> List[HardwareInfo]:
        return [self.get_mem_usage(), *self.get_cur_proc_mem(self.current_process)]
