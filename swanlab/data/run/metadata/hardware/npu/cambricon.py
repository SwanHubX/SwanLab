"""
@author: zeyi-lin
@file: cambrian.py
@time: 2025/03/27 12:32
@description: 寒武纪NPU信息采集
"""

import math
import platform
import subprocess

from ..type import HardwareFuncResult, HardwareInfoList, HardwareConfig, HardwareCollector as H
from ..utils import generate_key, random_index


def get_cambrian_npu_info() -> HardwareFuncResult:
    """
    获取寒武纪NPU信息，包括驱动版本、设备信息等
    """
    # cambrian芯片只支持Linux系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "npu": None}
    collector = None
    try:
        npu_map = map_npu()
        info["driver"] = npu_map[list(npu_map.keys())[0]]["driver"]
        # 将npu_map中的driver信息去掉
        for npu_id in npu_map:
            npu_map[npu_id].pop("driver")
        info["npu"] = npu_map
        collector = CambrianCollector(npu_map)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_npu() -> dict:
    """
    列出所有NPU设备，并包含芯片的映射关系
    返回值样例：
    {"0": { "name": "Cambrian-910", "memory": 16, "driver": "1.1.0"}, "1": { "name": "Cambrian-910", "memory": 16, "driver": "1.1.0"}, ...}
    """
    output = subprocess.run(["cnmon", "info", "-m"], capture_output=True, check=True, text=True).stdout
    npu_map = {}
    lines = output.split("\n")[1:]
    
    for line in lines:
        line = line.split()
        try:
            if line[0] == "Card":
                npu_id = line[-1]
                if npu_id not in npu_map:
                    npu_map[npu_id] = {}
            if line[0] == "Product":
                npu_name = line[-1]
                npu_map[npu_id]["name"] = npu_name
            if line[0] == "Driver" and len(npu_map[npu_id]) == 1:
                    driver = line[-1]
                    npu_map[npu_id]["driver"] = driver
            if line[0] == "Total" and len(npu_map[npu_id]) == 2:
                    memory = line[-2]
                    npu_map[npu_id]["memory"] = int(memory) // 1024  # 单位为GB
        except Exception as e:
            continue

    return npu_map


class CambrianCollector(H):
    def __init__(self, npu_map):
        super().__init__()
        self.npu_map = npu_map
        # NPU Utilization (%)
        self.util_key = generate_key("npu.{npu_index}.ptc")
        util_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="NPU Utilization (%)",
        )
        self.per_util_configs = {}
        # NPU Memory Allocated (%)
        self.hbm_rate_key = generate_key("npu.{npu_index}.mem.ptc")
        hbm_rate_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="NPU Memory Allocated (%)",
        )
        self.per_hbm_configs = {}
        self.per_temp_configs = {}

        for npu_id in npu_map:
            metric_name = f"NPU {npu_id}"
            self.per_util_configs[metric_name] = util_config.clone(metric_name=metric_name)
            self.per_hbm_configs[metric_name] = hbm_rate_config.clone(metric_name=metric_name)

    def collect(self) -> HardwareInfoList:
        result: HardwareInfoList = []
        for key, value in self.get_usage().items():
            result.append(value)
        for key, value in self.get_hbm_usage().items():
            result.append(value)
        return result

    def get_usage(self) -> HardwareInfoList:
        """
        获取指定NPU设备的利用率
        """
        output = subprocess.run(
            ["cnmon", "info", "-u"],
            capture_output=True,
            text=True,
        ).stdout
        # 获取信息
        
        util_infos = {}
        npu_ids = []
        
        # 获取所有NPU ID
        for npu_id in self.npu_map:
            npu_ids.append(npu_id)
        
        index = 0
        lines = output.split("\n")
        for line in lines:
            if "mlu average" in line.lower():             
                util_infos[npu_ids[index]] = {
                    "key": self.util_key.format(npu_index=npu_ids[index]),
                    "name": f"NPU {npu_ids[index]} Utilization (%)",
                    "value": math.nan,
                    "config": self.per_util_configs[f"NPU {npu_ids[index]}"],
                }
                line = line.split(":")
                # 利用率的值在最后一个
                util = line[-1].replace("%", "").strip()
                if util.isdigit():
                    util_infos[npu_ids[index]]['value'] = float(util)
                index += 1
                continue
                
        return util_infos
    
    def get_hbm_usage(self) -> HardwareInfoList:
        """
        获取指定NPU设备的HBM利用率
        """
        output = subprocess.run(
            ["cnmon", "info", "-m"],
            capture_output=True,
            text=True,
        ).stdout
        # 获取信息
        
        util_infos = {}
        npu_ids = []
        
        # 获取所有NPU ID
        for npu_id in self.npu_map:
            npu_ids.append(npu_id)
    
        index = 0
        lines = output.split("\n")
        for line_index, line in enumerate(lines):
            # 如果包含mlu average，则表示该行是NPU的利用率
            if "physical memory usage" in line.lower():
                if "used" in lines[line_index + 2].lower():
                    used_line = lines[line_index + 2]
                    # 初始化NPU的利用率
                    util_infos[npu_ids[index]] = {
                        "key": self.hbm_rate_key.format(npu_index=npu_ids[index]),
                        "name": f"NPU {npu_ids[index]} Memory Allocated (%)",
                        "value": math.nan,
                        "config": self.per_hbm_configs[f"NPU {npu_ids[index]}"],
                    }
                    used_line = used_line.split(":")
                    # 利用率的值在最后一个
                    util = used_line[-1].replace("MiB", "").strip()
                    if util.isdigit():
                        util_infos[npu_ids[index]]['value'] = float(util)
                    index += 1
                    continue
                
        return util_infos