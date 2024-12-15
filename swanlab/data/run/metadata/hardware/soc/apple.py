"""
@author: cunyue
@file: apple.py
@time: 2024/12/3 20:10
@description: 苹果芯片信息采集
"""

import json
import multiprocessing
import subprocess

import psutil
from swankit.env import is_macos

from ..type import HardwareFuncResult, HardwareInfoList, HardwareCollector as H
from ..utils import CpuCollector as C, MemoryCollector as M


def get_apple_chip_info() -> HardwareFuncResult:
    if not is_macos():
        return None, None
    info = {"cpu": None, "gpu": None, "memory": None, "type": None}
    # 使用system_profiler命令以JSON格式获取GPU信息
    try:
        result = subprocess.run(["system_profiler", "SPHardwareDataType", "-json"], capture_output=True, text=True)
        hardware_info = json.loads(result.stdout)["SPHardwareDataType"][0]
        chip_type = hardware_info.get("chip_type", None)
        # 早期intel芯片的机器
        if chip_type is None:
            chip_type = hardware_info.get("cpu_type", None)
        if chip_type is None:
            raise Exception("Can't get apple chip name")
        memory = hardware_info["physical_memory"]
        memory = str(memory).lower().replace("gb", "")
        # TODO: 获取GPU信息
        info["type"] = chip_type
        info["memory"] = memory
    except Exception:  # noqa
        return None, None
    try:
        info["cpu"] = multiprocessing.cpu_count()
    except Exception:  # noqa
        pass
    return info, AppleChipCollector()


class AppleChipCollector(H, C, M):
    def __init__(self):
        super().__init__()
        self.current_process = psutil.Process()

    def collect(self) -> HardwareInfoList:
        return [
            self.get_cpu_usage(),
            # *self.get_per_cpu_usage(),
            self.get_cur_proc_thds_num(self.current_process),
            self.get_mem_usage(),
            *self.get_cur_proc_mem(self.current_process),
        ]
