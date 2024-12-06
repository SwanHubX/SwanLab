"""
@author: cunyue
@file: memory.py
@time: 2024/12/3 20:13
@description: 内存信息采集
"""

import psutil

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult


def get_memory_size() -> HardwareFuncResult:
    """获取内存大小"""
    try:
        # 获取系统总内存大小
        mem = psutil.virtual_memory()
        total_memory = round(mem.total / (1024**3))  # 单位为GB
        return total_memory, []
    except Exception:  # noqa
        return None, []
