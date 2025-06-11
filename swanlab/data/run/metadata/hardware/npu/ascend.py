"""
@author: cunyue
@file: ascend.py
@time: 2024/11/26 12:32
@description: 华为昇腾NPU信息采集
"""

import math
import os
import platform
import subprocess
from typing import Optional

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfo, HardwareInfoList
from ..utils import generate_key, random_index


def get_ascend_npu_info() -> HardwareFuncResult:
    """
    获取华为昇腾NPU信息，包括驱动版本、设备信息等
    目前的信息统计粒度只到NPU ID级别，没有到Chip ID级别
    """
    # ascend芯片只支持Linux系统
    if platform.system() != "Linux":
        return None, None
    # /dev目录下没有davinci*设备文件，跳过
    # 其实理论上davinci后接数字，代表此设备id，但是官方文档也没明确写，以防万一还是不这么干了
    if not list(filter(lambda x: x.startswith("davinci"), os.listdir("/dev"))):
        return None, None
    info = {"driver": None, "npu": None, "cann": None}
    collector = None
    try:
        # 获取NPU驱动版本
        info["driver"] = get_version()
        # 获取CANN版本
        info["cann"] = get_cann_version()

        # 获取所有NPU设备ID
        npu_map = map_npu()
        hbm_value: int = 0
        for npu_id in npu_map:
            for chip_id in npu_map[npu_id]:
                chip_info = npu_map[npu_id][chip_id]
                usage = get_chip_usage(npu_id, chip_id)
                if info["npu"] is None:
                    info["npu"] = {}
                if npu_id not in info["npu"]:
                    info["npu"][npu_id] = {}
                info["npu"][npu_id][chip_id] = {**chip_info, "usage": usage}
                hbm_value = int(usage.get("hbm", 0)) if usage else 0
        hbm_value = hbm_value * 1024  # 转换为MB
        collector = AscendCollector(npu_map, hbm_value)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def get_version() -> str:
    result = subprocess.run(["npu-smi", "-v"], capture_output=True, check=True, text=True)
    return result.stdout.split(":")[-1].strip()


def get_cann_version() -> Optional[str]:
    """
    从 ascend_toolkit_install.info 文件中提取 Ascend-cann-toolkit 的版本号
    如果文件不存在或者格式不正确，返回 None
    """
    try:
        arch = platform.machine()
        if arch not in ("aarch64", "x86_64"):
            return None

        path = f"/usr/local/Ascend/ascend-toolkit/latest/{arch}-linux/ascend_toolkit_install.info"

        # 使用grep直接提取版本号，避免读取整个文件
        cmd = f"grep -m 1 '^version=' {path} | cut -d'=' -f2"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        return proc.stdout.strip() if proc.returncode == 0 else None

    except Exception:  # noqa
        return None


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
                usage["hbm"] = str(round(int(hbm) / 1024))  # 单位为GB
            break
    return usage


class AscendCollector(H):
    def __init__(self, npu_map, max_hbm_value: int):
        super().__init__()
        self.npu_map = npu_map
        self.max_hbm_value = max_hbm_value
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
        # NPU Memory Allocated (MB)
        self.hbm_value_key = generate_key("npu.{npu_index}.mem.value")
        hbm_value_config = HardwareConfig(
            y_range=(0, max_hbm_value),
            chart_index=random_index(),
            chart_name="NPU Memory Allocated (MB)",
        )
        self.per_hbm_value_configs = {}

        # NPU Temperature (℃)
        self.temp_key = generate_key("npu.{npu_index}.temp")
        temp_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="NPU Temperature (℃)",
        )
        self.per_temp_configs = {}
        # NPU Power Usage (W)
        self.power_key = generate_key("npu.{npu_index}.power")
        power_config = HardwareConfig(
            y_range=(0, None),
            chart_index=random_index(),
            chart_name="NPU Power Usage (W)",
        )
        self.per_power_config = {}

        for npu_id in npu_map:
            for chip_id in npu_map[npu_id]:
                metric_name = f"NPU {npu_id}-{chip_id}"
                self.per_util_configs[metric_name] = util_config.clone(metric_name=metric_name)
                self.per_hbm_configs[metric_name] = hbm_rate_config.clone(metric_name=metric_name)
                self.per_hbm_value_configs[metric_name] = hbm_value_config.clone(metric_name=metric_name)
                self.per_temp_configs[metric_name] = temp_config.clone(metric_name=metric_name)
                self.per_power_config[metric_name] = power_config.clone(metric_name=metric_name)

    def collect(self) -> HardwareInfoList:
        result: HardwareInfoList = []
        for npu_id in self.npu_map:
            for chip_id in self.npu_map[npu_id]:
                result.extend(self.get_usage(npu_id, chip_id))
                result.append(self.get_chip_temp(npu_id, chip_id))
                result.append(self.get_chip_power(npu_id, chip_id))
        return result

    def get_usage(self, npu_id: str, chip_id: str) -> HardwareInfoList:
        """
        获取指定NPU设备的芯片HBM的用量信息和利用率
        """
        output = subprocess.run(
            ["npu-smi", "info", "-t", "usages", "-i", npu_id, "-c", chip_id],
            capture_output=True,
            text=True,
        ).stdout
        # 格式化获取NPU ID和芯片ID
        _id, metric_name = self.get_label(npu_id, chip_id)
        # 获取信息
        util_info = {
            "key": self.util_key.format(npu_index=_id),
            "name": f"{metric_name} Utilization (%)",
            "value": math.nan,
            "config": self.per_util_configs[metric_name],
        }
        hbm_info = {
            "key": self.hbm_rate_key.format(npu_index=_id),
            "name": f"{metric_name} Memory Allocated (%)",
            "value": math.nan,
            "config": self.per_hbm_configs[metric_name],
        }
        hbm_value_info = {
            "key": self.hbm_value_key.format(npu_index=_id),
            "name": f"{metric_name} Memory Allocated (MB)",
            "value": math.nan,
            "config": self.per_hbm_value_configs[metric_name],
        }
        for line in output.split("\n"):
            if "aicore usage rate" in line.lower():
                line = line.split(":")
                # 利用率的值在最后一个
                util = line[-1].strip()
                if util.isdigit():
                    util_info["value"] = float(util)
                continue

            if "hbm usage rate" in line.lower():
                line = line.split(":")
                # HBM Capacity的值在最后一个
                hbm = line[-1].strip()
                if hbm.isdigit():
                    hbm_info["value"] = float(hbm)
                    hbm_value_info["value"] = float(hbm) * self.max_hbm_value * 0.01
                continue
        return [util_info, hbm_info, hbm_value_info]

    def get_chip_temp(self, npu_id: str, chip_id: str) -> HardwareInfo:
        """
        获取芯片温度
        """
        output = subprocess.run(
            ["npu-smi", "info", "-t", "temp", "-i", npu_id, "-c", chip_id],
            capture_output=True,
            text=True,
        ).stdout.strip()
        temp = float(output.split(":")[-1].strip())
        _id, metric_name = self.get_label(npu_id, chip_id)
        return {
            "key": self.temp_key.format(npu_index=_id),
            "name": f"{metric_name} Temperature (℃)",
            "value": temp,
            "config": self.per_temp_configs[metric_name],
        }

    def get_chip_power(self, npu_id: str, chip_id: str) -> HardwareInfo:
        """
        获取芯片功耗信息
        """
        output = subprocess.run(
            ["npu-smi", "info", "-t", "power", "-i", npu_id, "-c", chip_id],
            capture_output=True,
            text=True,
        ).stdout.strip()
        power = float(output.split(":")[-1].strip())
        _id, metric_name = self.get_label(npu_id, chip_id)
        return {
            "key": self.power_key.format(npu_index=_id),
            "name": f"{metric_name} Power Usage (W)",
            "value": power,
            "config": self.per_power_config[metric_name],
        }

    @staticmethod
    def get_label(npu_id: str, chip_id: str):
        _id = f"{npu_id}-{chip_id}"
        return _id, f"NPU {_id}"
