"""
@author: cunyue
@file: apple.py
@time: 2026/3/31 01:51
@description: 苹果统一芯片（Apple Silicon）相关的硬件供应商信息和功能实现
对于原本的Intel架构的Mac，我们在cpu.py和memory.py中统一采集信息
"""

import json
import multiprocessing
import subprocess
import sys
from typing import Optional

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.typings.run.system import AppleSiliconSnapshot
from swanlab.sdk.typings.run.system.hardware_vendor import AppleSiliconProtocol
from swanlab.sdk.utils.helper import catch_and_return_none


class Apple(AppleSiliconProtocol):
    """
    苹果统一芯片（Apple Silicon）相关的硬件供应商信息和功能实现
    """

    @staticmethod
    @catch_and_return_none(on_error=lambda e: console.debug("Failed to get Apple Silicon info: {}", e))
    def get() -> Optional[AppleSiliconSnapshot]:
        if sys.platform != "darwin":
            return None
        result = subprocess.run(
            ["system_profiler", "SPHardwareDataType", "-json"], capture_output=True, text=True, timeout=5
        )
        hardware_info = json.loads(result.stdout)["SPHardwareDataType"][0]
        chip_type = hardware_info.get("chip_type") or hardware_info.get("cpu_type")
        if chip_type is None:
            return None
        memory = str(hardware_info["physical_memory"]).lower().replace("gb", "").strip()
        cpu_count = multiprocessing.cpu_count()
        return AppleSiliconSnapshot(name=chip_type, memory=int(memory), memory_unit="GB", cpu_count=cpu_count)
