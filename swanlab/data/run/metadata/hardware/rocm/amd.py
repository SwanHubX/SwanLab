"""
@author: KMnO4-zx
@file: amd.py
@time: 2025-12-24
@description: AMD ROCm GPU 信息采集 (适配 SwanLab)
"""

import json
import math
import platform
import subprocess
import re
from typing import Optional, Tuple, List

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


# 用于标记 rocm-smi 是否可用（可能在某些 GPU 环境下会崩溃）
_ROCM_SMI_WORKING = True


def _mark_rocm_smi_failed():
    """标记 rocm-smi 不可用"""
    global _ROCM_SMI_WORKING
    _ROCM_SMI_WORKING = False


def get_amd_gpu_info() -> HardwareFuncResult:
    """
    获取 AMD GPU 信息，包括驱动版本、设备信息等
    支持 Linux (rocm-smi) 和 Windows (hipinfo)
    """
    system = platform.system()
    # Linux 使用 rocm-smi，Windows 使用 hipinfo
    if system not in ("Linux", "Windows"):
        return None, None

    info = {"driver": None, "gpu": None}
    collector = None
    try:
        driver, gpu_map = map_amd_gpu(system)
        info["driver"] = driver
        info["gpu"] = gpu_map
        max_mem_value = 0
        for gpu_id in gpu_map:
            # 提取显存数值，格式如 "64GB" -> 64
            try:
                mem_str = gpu_map[gpu_id]["memory"]
                # 移除 'GB' 并转为整数
                mem_value = int(float(re.sub(r"[^0-9.]", "", mem_str)))
                max_mem_value = max(max_mem_value, mem_value)
            except (ValueError, TypeError):
                continue

        # 将 GB 转为 MB 用于图表配置
        max_mem_value *= 1024
        collector = AMDCollector(gpu_map, max_mem_value, system)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_amd_gpu(system: str) -> Tuple[Optional[str], dict]:
    """
    获取 AMD GPU 信息，包括 ROCm 版本、设备列表
    返回示例:
    ROCm Version: 6.4.3
    {'0': {'name': 'AMD Radeon Graphics', 'memory': '64GB'}, ...}
    """
    if system == "Windows":
        return _map_amd_gpu_windows()

    # Linux 系统 - 使用 rocm-smi
    driver_version = "Unknown"
    gpu_map = _map_amd_gpu_via_rocm_smi()

    return driver_version, gpu_map


def _map_amd_gpu_windows() -> Tuple[Optional[str], dict]:
    """使用 hipinfo 获取 Windows 下的 AMD GPU 信息"""
    driver_version = "Unknown"
    gpu_map = {}
    try:
        result = subprocess.run(
            ["hipinfo"],
            capture_output=True,
            check=True,
            text=True,
        )
        output = result.stdout

        # 解析 hipinfo 输出
        # 查找所有设备
        device_blocks = re.split(r"^-{70,}", output, flags=re.MULTILINE)

        for block in device_blocks:
            # 解析设备 ID
            device_match = re.search(r"device#\s+(\d+)", block)
            if not device_match:
                continue
            device_id = device_match.group(1)

            # 解析设备名称
            name_match = re.search(r"Name:\s+(.+)", block)
            name = name_match.group(1).strip() if name_match else "AMD GPU"

            # 解析总显存 (格式: "107.87 GB")
            mem_match = re.search(r"totalGlobalMem:\s+([\d.]+)\s+GB", block)
            if mem_match:
                mem_gb = round(float(mem_match.group(1)))
                mem_str = f"{mem_gb}GB"
            else:
                mem_str = "0GB"

            # 解析 gcnArchName 获取架构信息
            arch_match = re.search(r"gcnArchName:\s+(\S+)", block)
            arch = arch_match.group(1) if arch_match else None

            gpu_map[device_id] = {"name": name, "memory": mem_str}
            if arch:
                gpu_map[device_id]["arch"] = arch

    except Exception:
        pass
    return driver_version, gpu_map


def _map_amd_gpu_via_rocm_smi() -> dict:
    """使用 rocm-smi 获取 GPU 设备信息"""
    global _ROCM_SMI_WORKING
    if not _ROCM_SMI_WORKING:
        return {}
    gpu_map = {}
    try:
        # 获取产品名称
        product_json = _run_rocm_smi_safe(["--showproductname", "--json"])
        # 获取显存信息 (VRAM)
        mem_json = _run_rocm_smi_safe(["--showmeminfo", "vram", "--json"])

        for card_key in product_json.keys():
            # card_key 通常是 "card0", "card1"
            idx = card_key.replace("card", "")

            # 获取名称
            name = product_json.get(card_key, {}).get("Card Series", "AMD GPU")

            # 获取显存大小 (通常单位是 Bytes)
            mem_bytes_str = mem_json.get(card_key, {}).get("VRAM Total Memory (B)", "0")
            try:
                mem_gb = round(int(mem_bytes_str) / (1024**3))
                mem_str = f"{mem_gb}GB"
            except (ValueError, TypeError):
                mem_str = "0GB"

            gpu_map[str(idx)] = {"name": name, "memory": mem_str}

    except Exception:
        pass
    return gpu_map


def _run_rocm_smi(args: list) -> dict:
    """
    辅助函数：运行 rocm-smi 并返回解析后的 JSON
    注意：如果 rocm-smi 崩溃，会标记为不可用
    """
    global _ROCM_SMI_WORKING
    if not _ROCM_SMI_WORKING:
        return {}
    return _run_rocm_smi_safe(args)


def _run_rocm_smi_safe(args: list) -> dict:
    """
    安全运行 rocm-smi，不修改全局状态
    """
    try:
        cmd = ["rocm-smi"] + args
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=5)
        return json.loads(result.stdout)
    except subprocess.TimeoutExpired:
        return {}
    except (subprocess.CalledProcessError, json.JSONDecodeError):
        return {}
    except Exception:
        # 其他异常（如崩溃），返回空字典
        return {}


class AMDCollector(H):
    def __init__(self, gpu_map, max_mem_value, system: str = "Linux"):
        super().__init__()
        self.gpu_map = gpu_map
        self.max_mem_value = max_mem_value
        self._system = system

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
            self.per_util_configs[metric_name] = util_config.clone(metric_name=metric_name)
            self.per_memory_configs[metric_name] = memory_config.clone(metric_name=metric_name)
            self.per_mem_value_configs[metric_name] = mem_value_config.clone(metric_name=metric_name)
            self.per_temp_configs[metric_name] = temp_config.clone(metric_name=metric_name)
            self.per_power_configs[metric_name] = power_config.clone(metric_name=metric_name)

    def collect(self) -> HardwareInfoList:
        result: HardwareInfoList = []
        if self._system == "Windows":
            # Windows 使用 hipinfo，仅支持内存监控
            usage_methods = [
                self.get_memory_usage,
                self.get_mem_value_usage,
            ]
        else:
            # Linux 使用 rocm-smi
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
        获取指定 GPU 设备的利用率
        rocm-smi: "GPU use (%)"
        """
        util_infos = {}
        output_json = _run_rocm_smi(["--showuse", "--json"])
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")

            util_infos[gpu_id] = {
                "key": self.util_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Utilization (%)",
                "value": math.nan,
                "config": self.per_util_configs.get(f"GPU {gpu_id}"),
            }

            val = info.get("GPU use (%)")
            if val is not None:
                try:
                    util_infos[gpu_id]["value"] = float(val)
                except (ValueError, TypeError):
                    pass
        return util_infos

    def get_memory_usage(self) -> dict:
        """
        获取指定 GPU 设备的内存占用率
        rocm-smi: "GPU Memory use (%)" 或 "VRAM use (%)"
        Windows: hipinfo memInfo.total 和 memInfo.free 计算
        """
        mem_infos = {}

        if self._system == "Windows":
            # Windows 使用 hipinfo
            devices = self._run_hipinfo()
            for device in devices:
                gpu_id = device.get("gpu", "")
                mem_infos[gpu_id] = {
                    "key": self.memory_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Memory Allocated (%)",
                    "value": math.nan,
                    "config": self.per_memory_configs.get(f"GPU {gpu_id}"),
                }

                total_mem = device.get("total_mem_gb")
                free_mem = device.get("free_mem_gb")
                if total_mem is not None and free_mem is not None and total_mem > 0:
                    used_mem = total_mem - free_mem
                    mem_infos[gpu_id]["value"] = (used_mem / total_mem) * 100
            return mem_infos

        # Linux 系统 - 使用 rocm-smi
        output_json = _run_rocm_smi(["--showmemuse", "--json"])
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")
            mem_infos[gpu_id] = {
                "key": self.memory_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Memory Allocated (%)",
                "value": math.nan,
                "config": self.per_memory_configs.get(f"GPU {gpu_id}"),
            }

            val = info.get("GPU Memory use (%)") or info.get("VRAM use (%)")
            if val is not None:
                try:
                    mem_infos[gpu_id]["value"] = float(val)
                except (ValueError, TypeError):
                    pass
        return mem_infos

    def get_mem_value_usage(self) -> dict:
        """
        获取指定 GPU 设备的内存使用量（MB）
        rocm-smi: 通过百分比计算
        Windows: hipinfo memInfo.total 和 memInfo.free 计算
        """
        mem_value_infos = {}

        if self._system == "Windows":
            # Windows 使用 hipinfo
            devices = self._run_hipinfo()
            for device in devices:
                gpu_id = device.get("gpu", "")
                mem_value_infos[gpu_id] = {
                    "key": self.mem_value_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Memory Allocated (MB)",
                    "value": math.nan,
                    "config": self.per_mem_value_configs.get(f"GPU {gpu_id}"),
                }

                total_mem = device.get("total_mem_gb")
                free_mem = device.get("free_mem_gb")
                if total_mem is not None and free_mem is not None:
                    used_mem_gb = total_mem - free_mem
                    mem_value_infos[gpu_id]["value"] = used_mem_gb * 1024  # 转为 MB
            return mem_value_infos

        # Linux 系统 - 使用 rocm-smi
        output_json = _run_rocm_smi(["--showmemuse", "--json"])
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")
            mem_value_infos[gpu_id] = {
                "key": self.mem_value_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Memory Allocated (MB)",
                "value": math.nan,
                "config": self.per_mem_value_configs.get(f"GPU {gpu_id}"),
            }

            val = info.get("GPU Memory use (%)") or info.get("VRAM use (%)")

            if val is not None:
                try:
                    current_card_mem_str = self.gpu_map.get(gpu_id, {}).get("memory", "0GB")
                    current_card_mem_mb = int(float(re.sub(r"[^0-9.]", "", current_card_mem_str))) * 1024

                    mem_value_infos[gpu_id]["value"] = float(val) * 0.01 * current_card_mem_mb
                except (ValueError, TypeError):
                    pass
        return mem_value_infos

    def get_temperature_usage(self) -> dict:
        """
        获取指定 GPU 设备的温度(°C)
        rocm-smi: "Temperature (Sensor junction) (C)" 或 "Temperature (Sensor edge) (C)"
        """
        temp_infos = {}
        output_json = _run_rocm_smi(["--showtemp", "--json"])
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")

            temp_infos[gpu_id] = {
                "key": self.temp_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Temperature (°C)",
                "value": math.nan,
                "config": self.per_temp_configs.get(f"GPU {gpu_id}"),
            }

            val = info.get("Temperature (Sensor junction) (C)") or info.get("Temperature (Sensor edge) (C)")
            if val is not None:
                try:
                    temp_infos[gpu_id]["value"] = float(val)
                except (ValueError, TypeError):
                    pass
        return temp_infos

    def get_power_usage(self) -> dict:
        """
        获取指定 GPU 设备的功耗(W)
        rocm-smi: "Average Graphics Package Power (W)" 或 "Socket Power (W)"
        """
        power_infos = {}
        output_json = _run_rocm_smi(["--showpower", "--json"])
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")

            power_infos[gpu_id] = {
                "key": self.power_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Power (W)",
                "value": math.nan,
                "config": self.per_power_configs.get(f"GPU {gpu_id}"),
            }

            val = info.get("Average Graphics Package Power (W)") or info.get("Socket Power (W)")
            if val is None:
                for k, v in info.items():
                    if "Power" in k and "(W)" in k:
                        val = v
                        break

            if val is not None:
                try:
                    power_infos[gpu_id]["value"] = float(val)
                except (ValueError, TypeError):
                    pass
        return power_infos

    def _run_hipinfo(self) -> List[dict]:
        """运行 hipinfo 并解析为设备列表"""
        devices = []
        try:
            result = subprocess.run(
                ["hipinfo"],
                capture_output=True,
                check=True,
                text=True,
            )
            output = result.stdout
            # 按分隔线分割设备块
            device_blocks = re.split(r"^-{70,}", output, flags=re.MULTILINE)

            for block in device_blocks:
                device_data = {}
                # 解析设备 ID
                device_match = re.search(r"device#\s+(\d+)", block)
                if not device_match:
                    continue
                device_data["gpu"] = device_match.group(1)

                # 解析内存信息 (格式: "memInfo.total: 107.87 GB", "memInfo.free: 107.72 GB (100%)")
                total_match = re.search(r"memInfo\.total:\s+([\d.]+)\s+GB", block)
                free_match = re.search(r"memInfo\.free:\s+([\d.]+)\s+GB", block)
                if total_match:
                    device_data["total_mem_gb"] = float(total_match.group(1))
                if free_match:
                    device_data["free_mem_gb"] = float(free_match.group(1))

                devices.append(device_data)
        except Exception:
            pass
        return devices