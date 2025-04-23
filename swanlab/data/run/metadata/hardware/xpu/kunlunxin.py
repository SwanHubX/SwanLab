"""
@author: zeyi-lin
@file: kunlunxin.py
@time: 2025/04/14 12:32
@description: 昆仑芯xpu信息采集
"""

import math
import platform
import subprocess
from typing import Tuple, Optional

from ..type import HardwareFuncResult, HardwareInfoList, HardwareConfig, HardwareCollector as H
from ..utils import generate_key, random_index


def get_kunlunxin_xpu_info() -> HardwareFuncResult:
    """
    获取昆仑芯xpu信息，包括驱动版本、设备信息等
    """
    # kunlunxin芯片只支持Linux系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "xpu": None}
    collector = None
    try:
        driver, xpu_map = map_xpu()
        info["driver"] = driver
        info["xpu"] = xpu_map
        collector = KunlunxinCollector(xpu_map)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def get_xpu_driver() -> Tuple[Optional[str], dict]:
    """
    获取xpu驱动版本
    """
    output = subprocess.run(["xpu-smi", "-q"], capture_output=True, check=True, text=True).stdout
    lines = output.split("\n")
    for line in lines:
        if "Driver Version" in line:
            return line.split(":")[-1].strip()
    return None


def map_xpu() -> Tuple[Optional[str], dict]:
    """
    列出所有xpu设备，并包含芯片的映射关系
    返回值样例:
        driver: "1.1.0"
        xpu_map: {"0": { "name": "kunlunxin-910", "memory": 16, }, "1": { "name": "kunlunxin-910", "memory": 16, }, ...}
    """
    output = subprocess.run(["xpu-smi", "-m"], capture_output=True, check=True, text=True).stdout
    
    xpu_map = {}
    driver = get_xpu_driver()
    lines = output.split("\n")

    for line in lines:
        # 将line按空格分为列表
        try:
            xpu_map[line.split()[1]] = {
                "name": (line.split()[21][1:] + " " + line.split()[22][:-1]).strip(),
                "memory": int(line.split()[18])//1024,
            }
        except Exception:  # noqa
            continue
        
    return driver, xpu_map


class KunlunxinCollector(H):
    def __init__(self, xpu_map):
        super().__init__()
        self.xpu_map = xpu_map

        # xpu Utilization (%)
        self.util_key = generate_key("xpu.{xpu_index}.ptc")
        util_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="XPU Utilization (%)",
        )
        self.per_util_configs = {}

        # xpu Memory Allocated (%)
        self.memory_key = generate_key("xpu.{xpu_index}.mem.ptc")
        memory_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="XPU Memory Allocated (%)",
        )
        self.per_memory_configs = {}

        # xpu Temperature (°C)
        self.temp_key = generate_key("xpu.{xpu_index}.temp")
        temp_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="XPU Temperature (°C)",
        )
        self.per_temp_configs = {}

        # xpu Power (W)
        self.power_key = generate_key("xpu.{xpu_index}.power")
        power_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="XPU Power (W)",
        )
        self.per_power_configs = {}

        for xpu_id in self.xpu_map:
            metric_name = f"XPU {xpu_id}"
            self.per_util_configs[metric_name] = util_config.clone(metric_name=metric_name)
            self.per_memory_configs[metric_name] = memory_config.clone(metric_name=metric_name)
            self.per_temp_configs[metric_name] = temp_config.clone(metric_name=metric_name)
            self.per_power_configs[metric_name] = power_config.clone(metric_name=metric_name)

    def collect(self) -> HardwareInfoList:
        result: HardwareInfoList = []
        usage_methods = [
            self.get_utilization_usage,
            self.get_memory_usage,
            self.get_temperature_usage,
            self.get_power_usage,
        ]

        for method in usage_methods:
            result.extend(method().values())
        return result

    def get_utilization_usage(self) -> dict:
        """
        获取指定xpu设备的利用率
        """
        output = subprocess.run(["xpu-smi", "-m"], capture_output=True, text=True).stdout

        util_infos = {}
        xpu_ids = []

        # 获取所有xpu ID
        for xpu_id in self.xpu_map:
            xpu_ids.append(xpu_id)

        index = 0
        lines = output.split("\n")
        for line in lines:
            try: 
                
                util = line.split()[19]
                if util.isdigit():
                    util_infos[xpu_ids[index]] = {
                        "key": self.util_key.format(xpu_index=xpu_ids[index]),
                        "name": f"XPU {xpu_ids[index]} Utilization (%)",
                        "value": math.nan,
                        "config": self.per_util_configs[f"XPU {xpu_ids[index]}"],
                    }
                    util_infos[xpu_ids[index]]['value'] = float(util)
                    index += 1
            except Exception:  # noqa
                continue

        return util_infos

    def get_memory_usage(self) -> dict:
        """
        获取指定xpu设备的显存占用率
        """
        output = subprocess.run(["xpu-smi", "-m"], capture_output=True, text=True).stdout

        memory_infos = {}
        xpu_ids = []

        # 获取所有xpu ID
        for xpu_id in self.xpu_map:
            xpu_ids.append(xpu_id)

        index = 0
        lines = output.split("\n")
        for line in lines:
            try:
                memory = line.split()[17]
                if memory.isdigit():
                    memory_infos[xpu_ids[index]] = {
                        "key": self.memory_key.format(xpu_index=xpu_ids[index]),
                        "name": f"XPU {xpu_ids[index]} Memory Allocated (%)",
                        "value": math.nan,
                        "config": self.per_memory_configs[f"XPU {xpu_ids[index]}"],
                    }
                    memory_infos[xpu_ids[index]]['value'] = float(memory) / (self.xpu_map[xpu_ids[index]]['memory'] * 1024) * 100
                    index += 1
            except Exception:  # noqa
                continue

        return memory_infos

    def get_temperature_usage(self) -> dict:
        """
        获取指定xpu设备的温度(°C)
        """
        output = subprocess.run(["xpu-smi", "-m"], capture_output=True, text=True).stdout
        temp_infos = {}
        xpu_ids = []

        # 获取所有xpu ID
        for xpu_id in self.xpu_map:
            xpu_ids.append(xpu_id)

        index = 0
        lines = output.split("\n")
        for line in lines:
            try:
                temp = line.split()[4]
                if temp.isdigit():
                    temp_infos[xpu_ids[index]] = {
                        "key": self.temp_key.format(xpu_index=xpu_ids[index]),
                        "name": f"XPU {xpu_ids[index]} Temperature (°C)",
                        "value": math.nan,
                        "config": self.per_temp_configs[f"XPU {xpu_ids[index]}"],
                    }
                    temp_infos[xpu_ids[index]]['value'] = float(temp)
                    index += 1
            except Exception:  # noqa
                continue

        return temp_infos

    def get_power_usage(self) -> dict:
        """
        获取指定xpu设备的功耗(W)
        """
        output = subprocess.run(["xpu-smi", "-m"], capture_output=True, text=True).stdout
        power_infos = {}
        xpu_ids = []

        # 获取所有xpu ID
        for xpu_id in self.xpu_map:
            xpu_ids.append(xpu_id)

        index = 0
        lines = output.split("\n")
        for line in lines:
            try:
                power = line.split()[8]
                if power.isdigit():
                    power_infos[xpu_ids[index]] = {
                        "key": self.power_key.format(xpu_index=xpu_ids[index]),
                        "name": f"XPU {xpu_ids[index]} Power (W)",
                        "value": math.nan,
                        "config": self.per_power_configs[f"XPU {xpu_ids[index]}"],
                    }
                    power_infos[xpu_ids[index]]['value'] = float(power)
                    index += 1
            except Exception:  # noqa
                continue

        return power_infos
