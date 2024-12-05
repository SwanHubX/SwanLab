"""
@author: cunyue
@file: ascend.py
@time: 2024/11/26 12:32
@description: 华为昇腾NPU信息采集
"""

import os
import platform
import subprocess

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult


def get_ascend_npu_info() -> HardwareFuncResult:
    """
    获取华为昇腾NPU信息，包括驱动版本、设备信息等
    目前的信息统计粒度只到NPU ID级别，没有到Chip ID级别
    """
    # ascend芯片只支持Linux系统
    if platform.system() != "Linux":
        return None, []
    # /dev目录下没有davinci*设备文件，跳过
    # 其实理论上davinci后接数字，代表此设备id，但是官方文档也没明确写，以防万一还是不这么干了
    if not list(filter(lambda x: x.startswith("davinci"), os.listdir("/dev"))):
        return None, []
    info = {"driver": None, "npu": None}
    try:
        # 获取NPU驱动版本
        info["driver"] = get_version()
        # 获取所有NPU设备ID
        npu_map = map_npu()
        for npu_id in npu_map:
            for chip_id in npu_map[npu_id]:
                chip_info = npu_map[npu_id][chip_id]
                usage = get_chip_usage(npu_id, chip_id)
                if info["npu"] is None:
                    info["npu"] = {}
                if npu_id not in info["npu"]:
                    info["npu"][npu_id] = {}
                info["npu"][npu_id][chip_id] = {**chip_info, "usage": usage}
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, []
    return info, []


def get_version() -> str:
    result = subprocess.run(["npu-smi", "-v"], capture_output=True, check=True, text=True)
    return result.stdout.split(":")[-1].strip()


def map_npu() -> dict:
    """
    列出所有NPU设备，并包含芯片的映射关系

    """
    output = subprocess.run(["npu-smi", "info", "-m"], capture_output=True, check=True, text=True).stdout
    # npu_id -> chip_id -> {"id": chip_logic_id, "name": chip_name}
    npu_map = {}
    # output的第一行是表头，从第二行开始是数据，分别是NPU ID、芯片ID、芯片逻辑ID、芯片名称
    # 会包含Mcu芯片，这些芯片不参与计算，所以不存在逻辑ID，需要滤除
    lines = output.split("\n")[1:]
    for line in lines:
        line = line.split()
        if len(line) < 4:
            continue
        # 前三个值分别为NPU ID、芯片ID、芯片逻辑ID，最后可能会有多个值，是芯片名称
        npu_id, chip_id, chip_logic_id, *chip_name = line
        # 如果chip_logic_id不是数字，说明是Mcu芯片，不参与计算
        if not chip_logic_id.isdigit():
            continue
        chip_name = " ".join(chip_name)
        if npu_id not in npu_map:
            npu_map[npu_id] = {}
        npu_map[npu_id][chip_id] = {"id": chip_logic_id, "name": chip_name}
    return npu_map


def get_chip_usage(npu_id: str, chip_id: str):
    """
    获取某个NPU设备的芯片信息
    不再需要获取chip的名称
    """
    output = subprocess.run(
        ["npu-smi", "info", "-t", "usages", "-i", npu_id, "-c", chip_id],
        capture_output=True,
        text=True,
    ).stdout
    usage = {}
    # 以空格分隔，找到HBM Capacity所在的行
    for line in output.split("\n"):
        if "hbm capacity" in line.lower():
            line = line.split(":")
            # HBM Capacity的值在最后一个
            hbm = line[-1].strip()
            if hbm.isdigit():
                usage["hbm"] = str(round(int(hbm) / 1024))
            break
    return usage
