"""
@author: cunyue
@file: __init__.py.py
@time: 2024/12/3 20:14
@description: 硬件信息采集
"""

from typing import Callable, List, Any, Optional

from .cpu import get_cpu_info
from .gpu.nvidia import get_nvidia_gpu_info
from .memory import get_memory_size
from .npu.ascend import get_ascend_npu_info
from .soc.apple import get_apple_chip_info
from .type import HardwareFuncResult, HardwareMonitorFunc, HardwareInfo

__all__ = ["get_hardware_info", "HardwareMonitorFunc", "HardwareInfo"]


def get_hardware_info() -> HardwareFuncResult:
    """
    采集硬件信息，包括CPU、GPU、内存、硬盘等
    """
    monitor_funcs = []
    m = dec_hardware_func(get_memory_size, monitor_funcs)
    c = dec_hardware_func(get_cpu_info, monitor_funcs)
    nvidia = dec_hardware_func(get_nvidia_gpu_info, monitor_funcs)
    ascend = dec_hardware_func(get_ascend_npu_info, monitor_funcs)
    apple = dec_hardware_func(get_apple_chip_info, monitor_funcs)

    info = {
        "memory": m,
        "cpu": c,
        "gpu": {
            "nvidia": nvidia,
        },
        "npu": {
            "ascend": ascend,
        },
        "soc": {
            "apple": apple,
        },
    }
    return filter_none(info, fallback={}), monitor_funcs


def dec_hardware_func(
    func: Callable[[], HardwareFuncResult],
    monitor_funcs: List[HardwareMonitorFunc],
) -> Optional[Any]:
    """
    装饰器，用于记录硬件信息采集函数
    """
    x, y = func()
    monitor_funcs.extend(y)
    return x


def filter_none(data, fallback=None):
    """
    过滤掉字典中值为None的键值对，只对字典有效
    """
    if isinstance(data, dict):
        data = {k: v for k, v in data.items() if v is not None}
        if all(v is None for v in data.values()):
            return fallback
    return data
