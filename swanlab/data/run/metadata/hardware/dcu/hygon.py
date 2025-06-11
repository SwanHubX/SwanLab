"""
@author: caddiesnew
@file: hygon.py
@time: 2025-06-04 12:08:51
@description: Hygon DCU 信息采集
"""

import json
import math
import platform
import subprocess
from typing import Optional, Tuple

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


def get_hygon_dcu_info() -> HardwareFuncResult:
    """
    获取 Hygon DCU信息，包括驱动版本、设备信息等
    """
    # Hygon DCU默认为Linux系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "dcu": None}
    collector = None
    try:
        driver, dcu_map = map_hygon_dcu()
        info["driver"] = driver
        info["dcu"] = dcu_map
        max_mem_value = 0
        for dcu_id in dcu_map:
            mem_value = int(dcu_map[dcu_id]["memory"][:-2])
            max_mem_value = max(max_mem_value, mem_value)
        max_mem_value *= 1024
        collector = DCUCollector(dcu_map, max_mem_value)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_hygon_dcu() -> Tuple[Optional[str], dict]:
    """
    获取 Hygon DCU 信息，包括驱动版本、设备信息等，例如：
    Driver Version: 6.2.28
    {'0': {'name': 'K500SM_AI', 'memory': '64GB'}, '1': {'name': 'K500SM_AI', 'memory': '64GB',...}}
    """
    driver_info = subprocess.run(["hy-smi", "--showdriverversion"], capture_output=True, check=True, text=True).stdout
    driver_version = None
    for line in driver_info.split("\n"):
        if "Driver Version" in line:
            driver_version = line.split(":")[-1]
            break

    dcu_product_map_str = subprocess.run(
        ["hy-smi", "--showproductname", "--json"],
        capture_output=True,
        check=True,
        text=True,
    ).stdout
    dcu_product_map = json.loads(dcu_product_map_str)

    dcu_available_mem_str = subprocess.run(
        ["hy-smi", "--showmemavailable", "--json"],
        capture_output=True,
        check=True,
        text=True,
    ).stdout
    dcu_available_mem_map = json.loads(dcu_available_mem_str)

    dcu_map = {}
    for idx, device_id in enumerate(dcu_product_map.keys()):
        product_info = dcu_product_map.get(device_id, {})
        device_name = product_info.get("Card Series", "Unknown")

        available_info = dcu_available_mem_map.get(device_id, {})
        available_mem = available_info.get("Available memory size (MiB)", "0")

        try:
            total_mem_gb = round(int(available_mem) / 1024)
            mem_str = f"{total_mem_gb}GB"
        except (ValueError, TypeError):
            mem_str = "0GB"

        dcu_map[str(idx)] = {"name": device_name, "memory": mem_str}
    return driver_version, dcu_map


class DCUCollector(H):
    def __init__(self, dcu_map, max_mem_value):
        super().__init__()
        self.dcu_map = dcu_map
        self.max_mem_value = max_mem_value

        # DCU Utilization (%)
        self.util_key = generate_key("dcu.{dcu_index}.pct")
        util_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="DCU Utilization (%)",
        )
        self.per_util_configs = {}

        # DCU Memory Allocated (%)
        self.memory_key = generate_key("dcu.{dcu_index}.mem.pct")
        memory_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="DCU Memory Allocated (%)",
        )
        self.per_memory_configs = {}

        # DCU Memory Allocated (MB)
        self.mem_value_key = generate_key("dcu.{dcu_index}.mem.value")
        mem_value_config = HardwareConfig(
            y_range=(0, self.max_mem_value),
            chart_index=random_index(),
            chart_name="DCU Memory Allocated (MB)",
        )
        self.per_mem_value_configs = {}

        # DCU Temperature (°C)
        self.temp_key = generate_key("dcu.{dcu_index}.temp")
        temp_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="DCU Temperature (°C)",
        )
        self.per_temp_configs = {}

        # DCU Power (W)
        self.power_key = generate_key("dcu.{dcu_index}.power")
        power_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="DCU Power (W)",
        )
        self.per_power_configs = {}

        for dcu_id in self.dcu_map:
            metric_name = f"DCU {dcu_id}"
            self.per_util_configs[metric_name] = util_config.clone(metric_name=metric_name)
            self.per_memory_configs[metric_name] = memory_config.clone(metric_name=metric_name)
            self.per_mem_value_configs[metric_name] = mem_value_config.clone(metric_name=metric_name)
            self.per_temp_configs[metric_name] = temp_config.clone(metric_name=metric_name)
            self.per_power_configs[metric_name] = power_config.clone(metric_name=metric_name)

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
        获取指定DCU设备的利用率
        """
        output_str = subprocess.run(
            ["hy-smi", "--showuse", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        util_infos = {}
        for idx, (dcu_key, dcu_info) in enumerate(output_json.items()):
            dcu_id = str(idx)

            util_infos[dcu_id] = {
                "key": self.util_key.format(dcu_index=dcu_id),
                "name": f"DCU {dcu_id} Utilization (%)",
                "value": math.nan,
                "config": self.per_util_configs[f"DCU {dcu_id}"],
            }
            dcu_utilization = dcu_info["DCU use (%)"]
            util_infos[dcu_id]["value"] = float(dcu_utilization)
        return util_infos

    def get_memory_usage(self) -> dict:
        """
        获取指定DCU设备的内存占用率
        """
        output_str = subprocess.run(
            ["hy-smi", "--showmemuse", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        mem_infos = {}

        for idx, (dcu_key, dcu_info) in enumerate(output_json.items()):
            dcu_id = str(idx)
            mem_infos[dcu_id] = {
                "key": self.memory_key.format(dcu_index=dcu_id),
                "name": f"DCU {dcu_id} Memory Allocated (%)",
                "value": math.nan,
                "config": self.per_memory_configs[f"DCU {dcu_id}"],
            }

            dcu_mem_use = dcu_info["DCU memory use (%)"]
            mem_infos[dcu_id]["value"] = float(dcu_mem_use)
        return mem_infos

    def get_mem_value_usage(self) -> dict:
        """
        获取指定DCU设备的内存使用量（MB）
        """
        output_str = subprocess.run(
            ["hy-smi", "--showmemuse", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        mem_value_infos = {}

        for idx, (dcu_key, dcu_info) in enumerate(output_json.items()):
            dcu_id = str(idx)
            mem_value_infos[dcu_id] = {
                "key": self.mem_value_key.format(dcu_index=dcu_id),
                "name": f"DCU {dcu_id} Memory Allocated (MB)",
                "value": math.nan,
                "config": self.per_mem_value_configs[f"DCU {dcu_id}"],
            }

            dcu_mem_use_rate = dcu_info["DCU memory use (%)"]
            mem_value_infos[dcu_id]["value"] = float(dcu_mem_use_rate) * 0.01 * self.max_mem_value
        return mem_value_infos

    def get_temperature_usage(self) -> dict:
        """
        获取指定DCU设备的温度(°C) (采集 Junction 核心温度)
        {
          "card0": {
            "Temperature (Sensor edge) (C)": "52.0",
            "Temperature (Sensor junction) (C)": "51.0", // collected
            "Temperature (Sensor mem) (C)": "52.0"
          }
        }
        """
        output_str = subprocess.run(
            ["hy-smi", "--showtemp", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        temp_infos = {}
        for idx, (dcu_key, dcu_info) in enumerate(output_json.items()):
            dcu_id = str(idx)

            temp_infos[dcu_id] = {
                "key": self.temp_key.format(dcu_index=dcu_id),
                "name": f"DCU {dcu_id} Temperature (°C)",
                "value": math.nan,
                "config": self.per_temp_configs[f"DCU {dcu_id}"],
            }
            dcu_temp = dcu_info["Temperature (Sensor junction) (C)"]
            temp_infos[dcu_id]["value"] = float(dcu_temp)
        return temp_infos

    def get_power_usage(self) -> dict:
        """
        获取指定DCUhy-smi --设备的功耗(W)
        """
        output_str = subprocess.run(
            ["hy-smi", "--showpower", "--json"],
            capture_output=True,
            text=True,
        ).stdout
        output_json = json.loads(output_str)
        power_infos = {}
        for idx, (dcu_key, dcu_info) in enumerate(output_json.items()):
            dcu_id = str(idx)

            power_infos[dcu_id] = {
                "key": self.power_key.format(dcu_index=dcu_id),
                "name": f"DCU {dcu_id} Power (W)",
                "value": math.nan,
                "config": self.per_power_configs[f"DCU {dcu_id}"],
            }
            dcu_power = dcu_info["Average Graphics Package Power (W)"]
            power_infos[dcu_id]["value"] = float(dcu_power)
        return power_infos
