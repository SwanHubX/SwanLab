"""
@author: zeyi-lin
@file: cambricon.py
@time: 2025/03/27 12:32
@description: 寒武纪mlu信息采集
"""

import math
import platform
import subprocess
from typing import Any, Dict, Optional, Tuple

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


def get_cambricon_mlu_info() -> HardwareFuncResult:
    """
    获取寒武纪mlu信息，包括驱动版本、设备信息等
    """
    # cambricon芯片只支持Linux系统
    if platform.system() != "Linux":
        return None, None

    info: Dict[str, Any] = {"driver": None, "mlu": None}
    collector = None
    try:
        driver, mlu_map = map_mlu()
        info["driver"] = driver
        info["mlu"] = mlu_map
        max_mem_value = 0
        for mlu_id in mlu_map:
            mlu_mem = int(mlu_map[mlu_id].get("memory", 0))
            max_mem_value = max(max_mem_value, mlu_mem)
        max_mem_value = max_mem_value * 1024  # 转换为MB
        collector = CambriconCollector(mlu_map, max_mem_value)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_mlu() -> Tuple[Optional[str], dict]:
    """
    列出所有mlu设备，并包含芯片的映射关系
    返回值样例:
        driver: "1.1.0"
        mlu_map: {"0": { "name": "cambricon-910", "memory": 16, }, "1": { "name": "cambricon-910", "memory": 16, }, ...}
    """
    output = subprocess.run(["cnmon", "info", "-m"], capture_output=True, check=True, text=True).stdout
    mlu_map = {}
    driver = None
    lines = output.split("\n")[1:]

    mlu_id = None
    for line in lines:
        line = line.split()
        try:
            if line[0] == "Card":
                mlu_id = line[-1]
                # 获取mlu的ID
                if mlu_id not in mlu_map:
                    mlu_map[mlu_id] = {}
            if mlu_id is None:
                continue
            # 获取mlu的名称
            if line[0] == "Product":
                # 如果没有mlu_id，则跳过（解析错误）
                mlu_name = line[-1]
                mlu_map[mlu_id]["name"] = mlu_name
            # 获取mlu的驱动版本
            if line[0] == "Driver" and len(mlu_map[mlu_id]) == 1 and driver is None:
                driver = line[-1]
            # 获取mlu的内存大小
            if line[0] == "Total" and len(mlu_map[mlu_id]) == 1:
                memory = line[-2]
                mlu_map[mlu_id]["memory"] = str(int(memory) // 1024)  # 单位为GB
        except Exception:  # noqa
            continue

    return driver, mlu_map


class CambriconCollector(H):
    def __init__(self, mlu_map, max_mem_value):
        super().__init__()
        self.mlu_map = mlu_map
        self.max_mem_value = max_mem_value

        # mlu Utilization (%)
        self.util_key = generate_key("mlu.{mlu_index}.ptc")
        util_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="MLU Utilization (%)",
        )
        self.per_util_configs = {}

        # mlu Memory Allocated (%)
        self.memory_key = generate_key("mlu.{mlu_index}.mem.ptc")
        memory_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="MLU Memory Allocated (%)",
        )
        self.per_memory_configs = {}

        # mlu Memory Allocated (MB)
        self.mem_value_key = generate_key("mlu.{mlu_index}.mem.value")
        mem_value_config = HardwareConfig(
            y_range=(0, self.max_mem_value),
            chart_index=random_index(),
            chart_name="MLU Memory Allocated (MB)",
        )
        self.per_mem_value_configs = {}

        # mlu Temperature (°C)
        self.temp_key = generate_key("mlu.{mlu_index}.temp")
        temp_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="MLU Temperature (°C)",
        )
        self.per_temp_configs = {}

        # mlu Power (W)
        self.power_key = generate_key("mlu.{mlu_index}.power")
        power_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="MLU Power (W)",
        )
        self.per_power_configs = {}

        for mlu_id in self.mlu_map:
            metric_name = f"MLU {mlu_id}"
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
        获取指定mlu设备的利用率
        """
        output = subprocess.run(["cnmon", "info", "-u"], capture_output=True, text=True).stdout

        util_infos = {}
        mlu_ids = []

        # 获取所有mlu ID
        for mlu_id in self.mlu_map:
            mlu_ids.append(mlu_id)

        index = 0
        lines = output.strip().split("\n")
        for line in lines:
            line = line.strip()
            if "mlu average" in line.lower():
                util_infos[mlu_ids[index]] = {
                    "key": self.util_key.format(mlu_index=mlu_ids[index]),
                    "name": f"MLU {mlu_ids[index]} Utilization (%)",
                    "value": math.nan,
                    "config": self.per_util_configs[f"MLU {mlu_ids[index]}"],
                }
                # 获得此mlu的利用率数值
                util = line.split(":")[-1].replace("%", "").strip()
                if util.isdigit():
                    util_infos[mlu_ids[index]]["value"] = float(util)
                index += 1
                continue
        return util_infos

    def get_memory_usage(self) -> dict:
        """
        获取指定mlu设备的显存占用率
        """
        output = subprocess.run(
            ["cnmon", "info", "-m"],
            capture_output=True,
            text=True,
        ).stdout

        memory_infos = {}
        mlu_ids = []

        # 获取所有mlu ID
        for mlu_id in self.mlu_map:
            mlu_ids.append(mlu_id)

        index = 0
        lines = output.strip().split("\n")
        for line_idx, line in enumerate(lines):
            # 如果包含mlu average，则表示该行是mlu的利用率
            if "physical memory usage" in line.lower():
                if "used" in lines[line_idx + 2].lower():
                    used_line = lines[line_idx + 2]
                    # 初始化mlu的利用率
                    memory_infos[mlu_ids[index]] = {
                        "key": self.memory_key.format(mlu_index=mlu_ids[index]),
                        "name": f"MLU {mlu_ids[index]} Memory Allocated (%)",
                        "value": math.nan,
                        "config": self.per_memory_configs[f"MLU {mlu_ids[index]}"],
                    }
                    # 获得此mlu的显存数值（MiB）
                    memory = used_line.split(":")[-1].replace("MiB", "").strip()
                    if memory.isdigit():
                        # 计算mlu显存占用率
                        memory_infos[mlu_ids[index]]["value"] = (
                            float(memory) / (float(self.mlu_map[mlu_ids[index]]["memory"]) * 1024) * 100
                        )
                    index += 1
                    continue
        return memory_infos

    def get_mem_value_usage(self) -> dict:
        """
        获取指定mlu设备的显存使用量 (MB)
        """
        output = subprocess.run(
            ["cnmon", "info", "-m"],
            capture_output=True,
            text=True,
        ).stdout

        mem_value_infos = {}
        mlu_ids = []

        # 获取所有mlu ID
        for mlu_id in self.mlu_map:
            mlu_ids.append(mlu_id)

        index = 0
        lines = output.strip().split("\n")
        for line_idx, line in enumerate(lines):
            # 如果包含mlu average，则表示该行是mlu的利用率
            if "physical memory usage" in line.lower():
                if "used" in lines[line_idx + 2].lower():
                    used_line = lines[line_idx + 2]
                    # 初始化mlu的利用率
                    mem_value_infos[mlu_ids[index]] = {
                        "key": self.mem_value_key.format(mlu_index=mlu_ids[index]),
                        "name": f"MLU {mlu_ids[index]} Memory Allocated (MB)",
                        "value": math.nan,
                        "config": self.per_mem_value_configs[f"MLU {mlu_ids[index]}"],
                    }
                    # 获得此mlu的显存数值（MiB）
                    memory = used_line.split(":")[-1].replace("MiB", "").strip()
                    if memory.isdigit():
                        mem_value_infos[mlu_ids[index]]["value"] = float(memory)
                    index += 1
                    continue
        return mem_value_infos

    def get_temperature_usage(self) -> dict:
        """
        获取指定mlu设备的温度(°C)
        """
        output = subprocess.run(["cnmon", "info", "-e"], capture_output=True, text=True).stdout
        temp_infos = {}
        mlu_ids = []

        # 获取所有mlu ID
        for mlu_id in self.mlu_map:
            mlu_ids.append(mlu_id)

        index = 0
        lines = output.strip().split("\n")
        for line_idx, line in enumerate(lines):
            if "temperature" in line.lower():
                if "chip" in lines[line_idx + 2].lower():
                    temp_line = lines[line_idx + 2]
                    temp_infos[mlu_ids[index]] = {
                        "key": self.temp_key.format(mlu_index=mlu_ids[index]),
                        "name": f"MLU {mlu_ids[index]} Temperature (°C)",
                        "value": math.nan,
                        "config": self.per_temp_configs[f"MLU {mlu_ids[index]}"],
                    }
                    # 获得此mlu的温度数值
                    temp = temp_line.split(":")[-1].replace("C", "").strip()
                    if temp.isdigit():
                        temp_infos[mlu_ids[index]]["value"] = float(temp)
                    index += 1
                    continue
        return temp_infos

    def get_power_usage(self) -> dict:
        """
        获取指定mlu设备的功耗(W)
        """
        output = subprocess.run(["cnmon", "info", "-p"], capture_output=True, text=True).stdout
        power_infos = {}
        mlu_ids = []

        # 获取所有mlu ID
        for mlu_id in self.mlu_map:
            mlu_ids.append(mlu_id)

        index = 0
        lines = output.strip().split("\n")
        for line_idx, line in enumerate(lines):
            if "power" in line.lower() and line_idx + 1 < len(lines):
                if "usage" in lines[line_idx + 1].lower():
                    power_line = lines[line_idx + 1]
                    power_infos[mlu_ids[index]] = {
                        "key": self.power_key.format(mlu_index=mlu_ids[index]),
                        "name": f"MLU {mlu_ids[index]} Power (W)",
                        "value": math.nan,
                        "config": self.per_power_configs[f"MLU {mlu_ids[index]}"],
                    }
                    # 获得此mlu的功耗数值
                    power = power_line.split(":")[-1].replace("W", "").strip()
                    if power.isdigit():
                        power_infos[mlu_ids[index]]["value"] = float(power)
                    index += 1
                    continue
        return power_infos
