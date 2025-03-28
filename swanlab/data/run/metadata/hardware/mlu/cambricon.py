"""
@author: zeyi-lin
@file: cambrian.py
@time: 2025/03/27 12:32
@description: 寒武纪mlu信息采集
"""

import math
import platform
import subprocess

from ..type import HardwareFuncResult, HardwareInfoList, HardwareConfig, HardwareCollector as H
from ..utils import generate_key, random_index


def get_cambrian_mlu_info() -> HardwareFuncResult:
    """
    获取寒武纪mlu信息，包括驱动版本、设备信息等
    """
    # cambrian芯片只支持Linux系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "mlu": None}
    collector = None
    try:
        mlu_map = map_mlu()
        info["driver"] = mlu_map[list(mlu_map.keys())[0]]["driver"]
        # 将mlu_map中的driver信息去掉
        for mlu_id in mlu_map:
            mlu_map[mlu_id].pop("driver")
        info["mlu"] = mlu_map
        collector = CambrianCollector(mlu_map)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_mlu() -> dict:
    """
    列出所有mlu设备，并包含芯片的映射关系
    返回值样例：
    {"0": { "name": "Cambrian-910", "memory": 16, "driver": "1.1.0"}, "1": { "name": "Cambrian-910", "memory": 16, "driver": "1.1.0"}, ...}
    """
    output = subprocess.run(["cnmon", "info", "-m"], capture_output=True, check=True, text=True).stdout
    mlu_map = {}
    lines = output.split("\n")[1:]
    
    for line in lines:
        line = line.split()
        try:
            if line[0] == "Card":
                mlu_id = line[-1]
                if mlu_id not in mlu_map:
                    mlu_map[mlu_id] = {}
            if line[0] == "Product":
                mlu_name = line[-1]
                mlu_map[mlu_id]["name"] = mlu_name
            if line[0] == "Driver" and len(mlu_map[mlu_id]) == 1:
                    driver = line[-1]
                    mlu_map[mlu_id]["driver"] = driver
            if line[0] == "Total" and len(mlu_map[mlu_id]) == 2:
                    memory = line[-2]
                    mlu_map[mlu_id]["memory"] = int(memory) // 1024  # 单位为GB
        except Exception as e:
            continue

    return mlu_map


class CambrianCollector(H):
    def __init__(self, mlu_map):
        super().__init__()
        self.mlu_map = mlu_map
        # mlu Utilization (%)
        self.util_key = generate_key("mlu.{mlu_index}.ptc")
        util_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="MLU Utilization (%)",
        )
        self.per_util_configs = {}
        # mlu Memory Allocated (%)
        self.hbm_rate_key = generate_key("mlu.{mlu_index}.mem.ptc")
        hbm_rate_config = HardwareConfig(
            y_range=(0, 100),
            chart_index=random_index(),
            chart_name="MLU Memory Allocated (%)",
        )
        self.per_hbm_configs = {}
        self.per_temp_configs = {}

        for mlu_id in mlu_map:
            metric_name = f"MLU {mlu_id}"
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
        获取指定mlu设备的利用率
        """
        output = subprocess.run(
            ["cnmon", "info", "-u"],
            capture_output=True,
            text=True,
        ).stdout
        # 获取信息
        
        util_infos = {}
        mlu_ids = []
        
        # 获取所有mlu ID
        for mlu_id in self.mlu_map:
            mlu_ids.append(mlu_id)
        
        index = 0
        lines = output.split("\n")
        for line in lines:
            if "mlu average" in line.lower():             
                util_infos[mlu_ids[index]] = {
                    "key": self.util_key.format(mlu_index=mlu_ids[index]),
                    "name": f"MLU {mlu_ids[index]} Utilization (%)",
                    "value": math.nan,
                    "config": self.per_util_configs[f"MLU {mlu_ids[index]}"],
                }
                line = line.split(":")
                # 利用率的值在最后一个
                util = line[-1].replace("%", "").strip()
                if util.isdigit():
                    util_infos[mlu_ids[index]]['value'] = float(util)
                index += 1
                continue
                
        return util_infos
    
    def get_hbm_usage(self) -> HardwareInfoList:
        """
        获取指定mlu设备的HBM利用率
        """
        output = subprocess.run(
            ["cnmon", "info", "-m"],
            capture_output=True,
            text=True,
        ).stdout
        # 获取信息
        
        util_infos = {}
        mlu_ids = []
        
        # 获取所有mlu ID
        for mlu_id in self.mlu_map:
            mlu_ids.append(mlu_id)
    
        index = 0
        lines = output.split("\n")
        for line_index, line in enumerate(lines):
            # 如果包含mlu average，则表示该行是mlu的利用率
            if "physical memory usage" in line.lower():
                if "used" in lines[line_index + 2].lower():
                    used_line = lines[line_index + 2]
                    # 初始化mlu的利用率
                    util_infos[mlu_ids[index]] = {
                        "key": self.hbm_rate_key.format(mlu_index=mlu_ids[index]),
                        "name": f"MLU {mlu_ids[index]} Memory Allocated (%)",
                        "value": math.nan,
                        "config": self.per_hbm_configs[f"MLU {mlu_ids[index]}"],
                    }
                    used_line = used_line.split(":")
                    # 利用率的值在最后一个
                    util = used_line[-1].replace("MiB", "").strip()
                    if util.isdigit():
                        util_infos[mlu_ids[index]]['value'] = float(util)
                    index += 1
                    continue
                
        return util_infos