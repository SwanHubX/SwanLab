"""
@author: zeyi-lin
@file: metax.py
@time: 2025-05-31 08:29:00
@description: Metax GPU信息采集
"""

import math
import platform
import subprocess
from typing import Optional, Tuple

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


def get_metax_gpu_info() -> HardwareFuncResult:
    """
    获取MetaX GPU信息，包括驱动版本、设备信息等
    """
    # MetaX GPU只支持Linux系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "maca": None, "gpu": None}
    collector = None
    try:
        driver, maca_version, gpu_map = map_metax_gpu()
        info["driver"] = driver
        info["maca"] = maca_version
        info["gpu"] = gpu_map
        max_mem_value: int = 0
        for gpu_id in gpu_map:
            mem_value = gpu_map[gpu_id]["memory"]
            if mem_value > max_mem_value:
                max_mem_value = mem_value
        max_mem_value = max_mem_value * 1024
        collector = MetaxCollector(gpu_map, max_mem_value)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_metax_gpu() -> Tuple[Optional[str], Optional[str], dict]:
    """
    获取 metax gpu 信息，包括驱动版本、设备信息等，例如：
    driver: '2.1.12'
    gpu_map: {"0": { "name": "MetaX C500”, "memory": "64GB"}, "1": { "name": "MetaX C500", "memory": "64GB"}, ...}
    """
    output_str = subprocess.run(["mx-smi"], capture_output=True, check=True, text=True).stdout

    driver = None
    maca_version = None
    gpu_map = {}
    index = 0

    output_str_line = output_str.split("\n")

    for line in output_str_line:
        if "mx-smi" in line:
            driver = line.split(" ")[-1]
        if "MACA" in line:
            maca_version = line.split(" ")[3]
        elif "MetaX" in line and "Management" not in line:
            gpu_map[index] = {}
            line_split = line.split(" ")
            line_split = [item for item in line_split if item != ""]
            gpu_map[index]["name"] = line_split[3]
        elif "MiB" in line:
            try:
                line_split = line.split(" ")
                line_split = [item for item in line_split if item != ""]
                gpu_memory = int(line_split[-4].split("/")[-1]) // 1024
                gpu_map[index]["memory"] = gpu_memory
                index += 1
            except Exception as e:
                pass

    return driver, maca_version, gpu_map


class MetaxCollector(H):
    def __init__(self, gpu_map, max_mem_value):
        super().__init__()
        self.gpu_map = gpu_map
        self.max_mem_value = max_mem_value

        # GPU Utilization (%)
        self.util_key = generate_key("gpu.{gpu_index}.pct")
        util_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="GPU Utilization (%)",
        )
        self.per_util_configs = {}

        # GPU Memory Allocated (%)
        self.memory_key = generate_key("gpu.{gpu_index}.mem.pct")
        memory_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="GPU Memory Allocated (%)",
        )
        self.per_memory_configs = {}

        # GPU Memory Allocated (MB)
        self.memory_value_key = generate_key("gpu.{gpu_index}.mem.value")
        memory_value_config = HardwareConfig(
            y_range=(0, self.max_mem_value),
            chart_index=random_index(),
            chart_name="GPU Memory Allocated (MB)",
        )
        self.per_memory_value_configs = {}

        # GPU Temperature (°C)
        self.temp_key = generate_key("gpu.{gpu_index}.temp")
        temp_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="GPU Temperature (°C)",
        )
        self.per_temp_configs = {}

        # GPU Power (W)
        self.power_key = generate_key("gpu.{gpu_index}.power")
        power_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="GPU Power (W)",
        )
        self.per_power_configs = {}

        for gpu_id in self.gpu_map:
            metric_name = f"GPU {gpu_id}"
            self.per_util_configs[metric_name] = util_config.clone(metric_name=metric_name)
            self.per_memory_configs[metric_name] = memory_config.clone(metric_name=metric_name)
            self.per_memory_value_configs[metric_name] = memory_value_config.clone(metric_name=metric_name)
            self.per_temp_configs[metric_name] = temp_config.clone(metric_name=metric_name)
            self.per_power_configs[metric_name] = power_config.clone(metric_name=metric_name)

    def collect(self) -> HardwareInfoList:
        result: HardwareInfoList = []
        usage_methods = [
            self.get_utilization_usage,
            self.get_memory_rate_usage,
            self.get_memory_value_usage,
            self.get_temperature_usage,
            self.get_power_usage,
        ]

        for method in usage_methods:
            result.extend(method().values())
        return result

    def get_utilization_usage(self) -> dict:
        """
        获取指定GPU设备的利用率
        """
        output_str = subprocess.run(["mx-smi", "--show-usage"], capture_output=True, check=True, text=True).stdout

        output_str_line = output_str.split("\n")
        usage_infos = {}
        index = 0

        for line in output_str_line:
            if "GPU" in line and "%" in line:
                gpu_id = index
                usage_line = line.split(" ")
                usage_line = [item for item in usage_line if item != ""]
                gpu_usage = float(usage_line[-2])

                usage_infos[gpu_id] = {
                    "key": self.util_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Utilization (%)",
                    "value": math.nan,
                    "config": self.per_util_configs[f"GPU {gpu_id}"],
                }

                usage_infos[gpu_id]["value"] = gpu_usage
                index += 1

        return usage_infos

    def get_memory_rate_usage(self) -> dict:
        """
        获取指定GPU设备的内存占用率
        """
        output_str = subprocess.run(
            ["mx-smi", "--show-memory"],
            capture_output=True,
            text=True,
        ).stdout
        output_str_line = output_str.split("\n")
        mem_infos = {}
        index = 0
        for line in output_str_line:
            if "vis_vram total" in line:
                gpu_id = index
                gpu_mem_total = line.split(" ")[-2]
            if "vis_vram used" in line:
                gpu_id = index
                gpu_mem_used = line.split(" ")[-2]

                mem_infos[gpu_id] = {
                    "key": self.memory_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Memory Allocated (%)",
                    "value": math.nan,
                    "config": self.per_memory_configs[f"GPU {gpu_id}"],
                }

                gpu_mem_usage_rate = float(gpu_mem_used) / float(gpu_mem_total) * 100
                mem_infos[gpu_id]["value"] = gpu_mem_usage_rate
                index += 1

        return mem_infos

    def get_memory_value_usage(self) -> dict:
        """
        获取指定GPU设备的内存占用值
        """
        output_str = subprocess.run(
            ["mx-smi", "--show-memory"],
            capture_output=True,
            text=True,
        ).stdout
        output_str_line = output_str.split("\n")
        mem_value_infos = {}
        index = 0
        for line in output_str_line:
            if "vis_vram used" in line:
                gpu_id = index
                gpu_mem_used = line.split(" ")[-2]
                mem_value_infos[gpu_id] = {
                    "key": self.memory_value_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Memory Allocated (MB)",
                    "value": math.nan,
                    "config": self.per_memory_value_configs[f"GPU {gpu_id}"],
                }
                gpu_mem_usage_value = float(gpu_mem_used) / 1024
                mem_value_infos[gpu_id]["value"] = gpu_mem_usage_value
                index += 1
        return mem_value_infos

    def get_temperature_usage(self) -> dict:
        """
        获取指定GPU设备的利用率
        """
        output_str = subprocess.run(["mx-smi", "--show-temperature"], capture_output=True, check=True, text=True).stdout

        output_str_line = output_str.split("\n")
        temp_infos = {}
        index = 0

        for line in output_str_line:
            if "hotspot" in line:
                gpu_id = index
                gpu_temp = line.split(" ")[-2]

                temp_infos[gpu_id] = {
                    "key": self.temp_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Temperature (°C)",
                    "value": math.nan,
                    "config": self.per_temp_configs[f"GPU {gpu_id}"],
                }
                temp_infos[gpu_id]["value"] = float(gpu_temp)
                index += 1

        return temp_infos

    def get_power_usage(self) -> dict:
        """
        获取指定GPU设备的功耗(W)
        """
        output_str = subprocess.run(["mx-smi", "--show-board-power"], capture_output=True, check=True, text=True).stdout

        output_str_line = output_str.split("\n")
        power_infos = {}
        index = 0

        for line in output_str_line:
            if "Total" in line and "W" in line:
                total_line = line.split(" ")
                total_line = [item for item in total_line if item != ""]
                total_power = float(total_line[-2])
                power_infos[index] = {
                    "key": self.power_key.format(gpu_index=index),
                    "name": f"GPU {index} Power (W)",
                    "value": math.nan,
                    "config": self.per_power_configs[f"GPU {index}"],
                }
                power_infos[index]["value"] = total_power
                index += 1

        return power_infos
