"""
@author: zeyi-lin
@file: cambrian.py
@time: 2025/03/27 12:32
@description: 寒武纪NPU信息采集
"""

import platform
import subprocess

from ..type import HardwareFuncResult


def get_cambrian_npu_info() -> HardwareFuncResult:
    """
    获取寒武纪NPU信息，包括驱动版本、设备信息等
    """
    # cambrian芯片只支持Linux系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "npu": None}
    try:
        npu_map = map_npu()
        info["driver"] = npu_map[list(npu_map.keys())[0]]["driver"]
        # 将npu_map中的driver信息去掉
        for npu_id in npu_map:
            npu_map[npu_id].pop("driver")
        info["npu"] = npu_map
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, None


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