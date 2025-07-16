"""
@author: cunyue
@file: __init__.py
@time: 2024/11/18 15:02
@description: 实验元信息采集
"""

from typing import Tuple, List

from swanlab.package import get_package_version
from swanlab.swanlab_settings import get_settings
from .conda import get_conda
from .cooperation import get_cooperation_info
from .hardware import *
from .requirements import get_requirements
from .runtime import get_runtime_info


def get_metadata(logdir: str = None) -> Tuple[dict, List[HardwareCollector]]:
    """
    采集实验的全部信息
    """
    settings = get_settings()
    # 1. 按照配置采集元信息
    hardware_info, monitor_funcs, runtime_info = {}, [], {}
    if settings.metadata_collect:
        # 1.1 采集软硬件信息
        if settings.collect_hardware:
            hardware_info, monitor_funcs = get_hardware_info()
        # 1.2 如果硬件监控被关闭，则monitor_funcs置空
        if not settings.hardware_monitor:
            monitor_funcs = []
        # 1.3 采集运行时信息
        if settings.collect_runtime:
            runtime_info = get_runtime_info()
    # 2. swanlab官方信息收集
    coop = get_cooperation_info()
    # 2.1 生成基本swanlab信息
    swanlab_info = {
        "version": get_package_version(),
        "_settings": settings.filter_changed_fields(),
        "_coop": coop,
        "_monitor": len(monitor_funcs),
    }
    if logdir is not None:
        swanlab_info["logdir"] = logdir
    # 2.2 如果_coop为空，则删除
    if not coop:
        del swanlab_info["_coop"]
    # 2.3 如果_settings为空，则删除
    if not swanlab_info["_settings"]:
        del swanlab_info["_settings"]
    # 3. 合并所有信息
    metadata = {
        **hardware_info,
        **runtime_info,
        "swanlab": swanlab_info,
    }
    return metadata, monitor_funcs


__all__ = [
    "get_metadata",
    "get_requirements",
    "get_conda",
    "get_cooperation_info",
    "HardwareInfo",
    "HardwareCollector",
]
