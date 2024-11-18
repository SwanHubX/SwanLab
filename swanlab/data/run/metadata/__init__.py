"""
@author: cunyue
@file: __init__.py
@time: 2024/11/18 15:02
@description: 实验元信息采集
"""

from swanlab.data.run.metadata.cooperation import get_cooperation_info
from swanlab.data.run.metadata.hardware import get_hardware_info
from swanlab.data.run.metadata.requirements import get_requirements
from swanlab.data.run.metadata.runtime import get_runtime_info


def get_metadata(logdir: str):
    """
    采集实验的全部信息
    """
    coop = get_cooperation_info()
    return {
        **get_hardware_info(),
        **get_runtime_info(),
        "swanlab": {
            "version": coop["swanlab"]["version"],
            "logdir": logdir,
            "_coop": coop,
        },
    }


__all__ = ["get_metadata", "get_requirements", "get_cooperation_info"]
