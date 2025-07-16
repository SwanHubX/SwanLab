"""
@author: caddiesnew
@file: moorethread.py
@time: 2025-05-20 16:39:17
@description: MooreThread GPU信息采集
"""

import json
import math
import platform
import subprocess
from typing import Optional, Tuple

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


def get_moorethreads_gpu_info() -> HardwareFuncResult:
    """
    获取MooreThreads GPU信息，包括驱动版本、设备信息等
    """
    # MooreThreads GPU只支持Linux系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "gpu": None}
    collector = None
    try:
        driver, gpu_map = map_moorethreads_gpu()
        info["driver"] = driver
        info["gpu"] = gpu_map
        max_mem_value = 0
        for gpu_id in gpu_map:
            mem_value = int(gpu_map[gpu_id]["memory"])
            max_mem_value = max(max_mem_value, mem_value)
        max_mem_value *= 1024
        collector = MTTCollector(gpu_map, max_mem_value)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_moorethreads_gpu() -> Tuple[Optional[str], dict]:
    """
    获取 Moore Threads GPU信息，包括驱动版本、设备信息等，例如：
    driver: '2.7.0'
    gpu_map: {"0": { "name": "MTT S4000”, "memory": "48"}, "1": { "name": "MTT S4000", "memory": "48"}, ...}
    """
    output_str = subprocess.run(
        ["mthreads-gmi", "-q", "--json"], capture_output=True, check=True, text=True
    ).stdout
    output_json = json.loads(output_str)
    driver = None
    gpu_map = {}

    if "Driver Version" in output_json.keys():
        driver = output_json["Driver Version"]

    for gpu_info in output_json["GPU"]:
        gpu_id = gpu_info["Index"]
        gpu_name = gpu_info["Product Name"]
        gpu_memory = gpu_info["FB Memory Usage"]["Total"]
        if gpu_memory.endswith("MiB"):
            total = int(gpu_memory[:-3]) // 1024
            gpu_memory = str(total)
        elif gpu_memory.endswith("GiB"):
            total = int(gpu_memory[:-3])
            gpu_memory = str(total)
        else:
            pass
        gpu_map[gpu_id] = {
            "name": gpu_name,
            "memory": gpu_memory,
        }
    return driver, gpu_map


class MTTCollector(H):
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
        self.mem_value_key = generate_key("gpu.{gpu_index}.mem.value")
        mem_value_config = HardwareConfig(
            y_range=(0, self.max_mem_value),
            chart_index=random_index(),
            chart_name="GPU Memory Allocated (MB)",
        )
        self.per_mem_value_configs = {}

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
            self.per_util_configs[metric_name] = util_config.clone(
                metric_name=metric_name
            )
            self.per_memory_configs[metric_name] = memory_config.clone(
                metric_name=metric_name
            )
            self.per_mem_value_configs[metric_name] = mem_value_config.clone(
                metric_name=metric_name
            )
            self.per_temp_configs[metric_name] = temp_config.clone(
                metric_name=metric_name
            )
            self.per_power_configs[metric_name] = power_config.clone(
                metric_name=metric_name
            )

    def collect(self) -> HardwareInfoList:
        result: HardwareInfoList = []
        usage_methods = [
            self.get_utilization_usage,
            self.get_memory_usage,
            self.get_mem_value_usage,
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
        output_str = subprocess.run(
            ["mthreads-gmi", "-q", "-d", "UTILIZATION", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        util_infos = {}
        for gpu_info in output_json["GPU"]:
            gpu_id = gpu_info["Index"]

            util_infos[gpu_id] = {
                "key": self.util_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Utilization (%)",
                "value": math.nan,
                "config": self.per_util_configs[f"GPU {gpu_id}"],
            }
            # 解析利用率
            gpu_utilization = gpu_info["Utilization"]["Gpu"]

            if isinstance(gpu_utilization, str) and gpu_utilization.endswith("%"):
                gpu_utilization = gpu_utilization[:-1]
            try:
                util_infos[gpu_id]["value"] = float(gpu_utilization)
            except ValueError:
                pass
        return util_infos

    def get_memory_usage(self) -> dict:
        """
        获取指定GPU设备的内存占用率
        """
        output_str = subprocess.run(
            ["mthreads-gmi", "-q", "-d", "MEMORY", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        mem_infos = {}
        for gpu_info in output_json["GPU"]:
            gpu_id = gpu_info["Index"]

            mem_infos[gpu_id] = {
                "key": self.memory_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Memory Allocated (%)",
                "value": math.nan,
                "config": self.per_memory_configs[f"GPU {gpu_id}"],
            }
            gpu_mem_info = gpu_info["FB Memory Usage"]
            gpu_mem_toal = gpu_mem_info["Total"].replace("MiB", "").strip()
            gpu_mem_used = gpu_mem_info["Used"].replace("MiB", "").strip()
            gpu_mem_usage_rate = float(gpu_mem_used) / float(gpu_mem_toal) * 100
            mem_infos[gpu_id]["value"] = gpu_mem_usage_rate
        return mem_infos

    def get_mem_value_usage(self) -> dict:
        """
        获取指定GPU设备的显存使用量
        """
        output_str = subprocess.run(
            ["mthreads-gmi", "-q", "-d", "MEMORY", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        mem_infos = {}
        for gpu_info in output_json["GPU"]:
            gpu_id = gpu_info["Index"]

            mem_infos[gpu_id] = {
                "key": self.mem_value_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Memory Allocated (MB)",
                "value": math.nan,
                "config": self.per_mem_value_configs[f"GPU {gpu_id}"],
            }
            gpu_mem_info = gpu_info["FB Memory Usage"]
            mem_infos[gpu_id]["value"] = float(
                gpu_mem_info["Used"].replace("MiB", "").strip()
            )
        return mem_infos

    def get_temperature_usage(self) -> dict:
        """
        获取指定GPU设备的温度(°C)
        """
        output_str = subprocess.run(
            ["mthreads-gmi", "-q", "-d", "TEMPERATURE", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        temp_infos = {}
        for gpu_info in output_json["GPU"]:
            gpu_id = gpu_info["Index"]

            temp_infos[gpu_id] = {
                "key": self.temp_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Temperature (°C)",
                "value": math.nan,
                "config": self.per_temp_configs[f"GPU {gpu_id}"],
            }
            # 解析温度
            gpu_temp = (
                gpu_info["Temperature"]["GPU Current Temp"].replace("C", "").strip()
            )
            temp_infos[gpu_id]["value"] = float(gpu_temp)
        return temp_infos

    def get_power_usage(self) -> dict:
        """
        获取指定GPU设备的功耗(W)
        """
        output_str = subprocess.run(
            ["mthreads-gmi", "-q", "-d", "POWER", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        power_infos = {}
        for gpu_info in output_json["GPU"]:
            gpu_id = gpu_info["Index"]

            power_infos[gpu_id] = {
                "key": self.power_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Power (W)",
                "value": math.nan,
                "config": self.per_power_configs[f"GPU {gpu_id}"],
            }
            gpu_power = (
                gpu_info["Power Readings"]["Power Draw "].strip().replace("W", "")
            )
            power_infos[gpu_id]["value"] = float(gpu_power)
        return power_infos
