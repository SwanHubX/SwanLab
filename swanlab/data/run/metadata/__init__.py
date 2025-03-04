"""
@author: cunyue
@file: __init__.py
@time: 2024/11/18 15:02
@description: 实验元信息采集
"""

from typing import Tuple, List

from swanlab.package import get_package_version
from .conda import get_conda
from .cooperation import get_cooperation_info
from .hardware import *
from .requirements import get_requirements
from .runtime import get_runtime_info


def get_metadata(logdir: str) -> Tuple[dict, List[HardwareCollector]]:
    """
    采集实验的全部信息
    """
    hardware_info, monitor_funcs = get_hardware_info()
    coop = get_cooperation_info()
    metadata = {
        **hardware_info,
        **get_runtime_info(),
        "swanlab": {
            "version": get_package_version(),
            "logdir": logdir,
            "_monitor": len(monitor_funcs),
            "_coop": coop,
        },
    }
    if not metadata['swanlab'].get("_coop"):
        metadata['swanlab'].pop("_coop")
    if not metadata['swanlab'].get("_monitor"):
        metadata['swanlab'].pop("_monitor")
    return metadata, monitor_funcs


__all__ = ["get_metadata", "get_requirements", "get_conda", "get_cooperation_info", "HardwareInfo", "HardwareCollector"]
