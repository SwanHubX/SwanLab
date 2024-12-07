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
from typing import Optional, List

import psutil
from swankit.callback.models import ColumnConfig

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult, HardwareInfo, HardwareCollector


def get_apple_chip_info() -> HardwareFuncResult:
    if "mac" not in platform.platform().lower():
        return None, None
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
        return None, None
    try:
        info["cpu"] = multiprocessing.cpu_count()
    except Exception:  # noqa
        pass
    return info, AppleChipCollector()


class AppleChipCollector(HardwareCollector):
    def __init__(self):
        self.cpu_config: ColumnConfig = {"y_range": (0, 100)}

    def collect(self) -> List[Optional[HardwareInfo]]:
        info = [self.get_cpu_usage()]
        return info

    def get_cpu_usage(self) -> HardwareInfo:
        value = psutil.cpu_percent(interval=1)
        return {
            "key": "apple_cpu_usage",
            "value": value,
            "name": "System CPU Utilization (%)",
            "config": self.cpu_config,
        }
