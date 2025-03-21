"""
@author: cunyue
@file: __init__.py.py
@time: 2024/12/3 20:14
@description: 硬件信息采集
"""

from typing import Callable, List, Any, Optional, Tuple

from .cpu import get_cpu_info
from .disk import get_disk_info
from .gpu.nvidia import get_nvidia_gpu_info
from .memory import get_memory_size
from .network import get_network_info
from .npu.ascend import get_ascend_npu_info
from .soc.apple import get_apple_chip_info
from .type import HardwareFuncResult, HardwareCollector, HardwareInfo

__all__ = ["get_hardware_info", "HardwareCollector", "HardwareInfo"]


def get_hardware_info() -> Tuple[Optional[Any], List[HardwareCollector]]:
    """
    采集硬件信息，包括CPU、GPU、内存、硬盘等
    """
    monitor_funcs = []
    # 我们希望计算芯片的信息放在最前面，前端展示用
    nvidia = dec_hardware_func(get_nvidia_gpu_info, monitor_funcs)
    ascend = dec_hardware_func(get_ascend_npu_info, monitor_funcs)
    apple = dec_hardware_func(get_apple_chip_info, monitor_funcs)
    c = dec_hardware_func(get_cpu_info, monitor_funcs)
    m = dec_hardware_func(get_memory_size, monitor_funcs)
    d = dec_hardware_func(get_disk_info, monitor_funcs)
    n = dec_hardware_func(get_network_info, monitor_funcs)

    info = {
        "memory": m,
        "cpu": c,
        "disk": d,
        "network": n,
        "gpu": {},
        "npu": {},
        "soc": {},
    }
    if nvidia is not None:
        info["gpu"]["nvidia"] = nvidia
    if ascend is not None:
        info["npu"]["ascend"] = ascend
    if apple is not None:
        info["soc"]["apple"] = apple

    return filter_none(info, fallback={}), monitor_funcs


def dec_hardware_func(
    func: Callable[[], HardwareFuncResult],
    monitor_funcs: List[HardwareCollector],
) -> Optional[Any]:
    """
    装饰器，用于记录硬件信息采集函数
    """
    x, y = func()
    if y:
        monitor_funcs.append(y)
    return x


def filter_none(data, fallback=None):
    """
    过滤掉字典中值为None的键值对，只对字典有效
    """
    if isinstance(data, dict):
        data = {k: v for k, v in data.items() if v is not None and v != {}}  # 过滤掉空字典
        if all(v is None for v in data.values()):
            return fallback
    return data
