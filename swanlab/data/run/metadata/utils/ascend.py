"""
@author: cunyue
@file: ascend.py
@time: 2024/11/26 12:32
@description: 华为昇腾NPU信息采集
"""

import subprocess


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
