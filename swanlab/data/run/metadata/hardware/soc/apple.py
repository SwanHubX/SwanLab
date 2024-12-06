"""
@author: cunyue
@file: apple.py
@time: 2024/12/3 20:10
@description: 苹果芯片信息采集
"""

import json
import multiprocessing
import platform
import subprocess

import psutil

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult, HardwareInfo
from swanlab.data.run.metadata.hardware.utils import hardware


def get_apple_chip_info() -> HardwareFuncResult:
    if "mac" not in platform.platform().lower():
        return None, []
    info = {"cpu": None, "gpu": None, "memory": None, "type": None}

    # 使用system_profiler命令以JSON格式获取GPU信息
    try:
        result = subprocess.run(["system_profiler", "SPHardwareDataType", "-json"], capture_output=True, text=True)
        gpu_name = json.loads(result.stdout)["SPHardwareDataType"][0]["chip_type"]
        memory = json.loads(result.stdout)["SPHardwareDataType"][0]["physical_memory"]
        memory = str(memory).lower().replace("gb", "")
        # TODO: 获取GPU信息
        info["type"] = gpu_name
        info["memory"] = memory
    except Exception:  # noqa
        return None, []
    try:
        info["cpu"] = multiprocessing.cpu_count()
    except Exception:  # noqa
        pass
    return info, [
        # get_cpu_usage
    ]


@hardware
def get_cpu_usage() -> HardwareInfo:
    usage = psutil.cpu_percent(interval=1)
    return {"key": "cpu_usage", "value": usage, "name": "System CPU Utilization (%)"}
