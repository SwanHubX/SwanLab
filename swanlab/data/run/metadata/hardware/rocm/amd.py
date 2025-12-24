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
import shutil
from typing import Optional, Tuple, List

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


# 检测可用的 SMI 工具: amd-smi (新版本) 或 rocm-smi (旧版本)
_AMD_SMI_AVAILABLE = shutil.which("amd-smi") is not None
_ROCM_SMI_AVAILABLE = shutil.which("rocm-smi") is not None

# 用于标记 amd-smi 是否可用 (可能在运行时发现不兼容)
_AMD_SMI_WORKING = True

# 使用 amd-smi 作为首选工具
_USE_AMD_SMI = _AMD_SMI_AVAILABLE


def _mark_amd_smi_failed():
    """标记 amd-smi 不可用，后续自动回退到 rocm-smi"""
    global _AMD_SMI_WORKING, _USE_AMD_SMI
    _AMD_SMI_WORKING = False
    _USE_AMD_SMI = False


def get_amd_gpu_info() -> HardwareFuncResult:
    """
    获取 AMD GPU 信息，包括驱动版本、设备信息等
    支持 Linux (amd-smi/rocm-smi) 和 Windows (hipinfo)
    """
    system = platform.system()
    # Linux 使用 amd-smi/rocm-smi，Windows 使用 hipinfo
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

    # Linux 系统
    # 1. 获取版本信息 (优先使用 amd-smi version，因为它包含详细的 ROCm 版本)
    driver_version = None

    if _USE_AMD_SMI:
        try:
            version_output = subprocess.run(
                ["amd-smi", "version"],
                capture_output=True,
                check=True,
                text=True,
                timeout=5,
            ).stdout
            # 解析 "ROCm version: 6.4.3"
            match = re.search(r"ROCm version:\s*([\d\.]+)", version_output)
            if match:
                driver_version = match.group(1)
            else:
                driver_version = "Unknown"
        except Exception:
            # amd-smi version 失败，标记为不可用
            _mark_amd_smi_failed()
            driver_version = "Unknown"

    # 2. 获取设备名称和显存信息
    gpu_map = {}

    if _USE_AMD_SMI:
        gpu_map = _map_amd_gpu_via_amd_smi()
        # 如果 amd-smi 获取失败且 rocm-smi 可用，尝试回退
        if not gpu_map and _ROCM_SMI_AVAILABLE:
            gpu_map = _map_amd_gpu_via_rocm_smi()
    elif _ROCM_SMI_AVAILABLE:
        gpu_map = _map_amd_gpu_via_rocm_smi()

    if not driver_version:
        driver_version = "Unknown"

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


def _map_amd_gpu_via_amd_smi() -> dict:
    """使用 amd-smi static 获取 GPU 设备信息"""
    global _AMD_SMI_WORKING
    if not _AMD_SMI_WORKING:
        return {}
    gpu_map = {}
    try:
        # 添加超时防止命令卡死
        result = subprocess.run(
            ["amd-smi", "static", "--json"],
            capture_output=True,
            check=True,
            text=True,
            timeout=5,
        )
        static_info = json.loads(result.stdout)

        for device in static_info:
            gpu_id = str(device.get("gpu", ""))
            # 获取设备名称
            asic_info = device.get("ASIC", {})
            name = asic_info.get("MARKET_NAME", "AMD GPU")

            # 获取显存大小 (VRAM SIZE 单位通常是 MB)
            vram_info = device.get("VRAM", {})
            vram_mb = vram_info.get("SIZE", 0)
            try:
                mem_gb = round(int(vram_mb) / 1024)
                mem_str = f"{mem_gb}GB"
            except (ValueError, TypeError):
                mem_str = "0GB"

            gpu_map[gpu_id] = {"name": name, "memory": mem_str}
    except subprocess.TimeoutExpired:
        # 命令超时，标记为失败并回退
        _mark_amd_smi_failed()
    except (subprocess.CalledProcessError, json.JSONDecodeError):
        # 命令执行失败，标记为失败并回退
        _mark_amd_smi_failed()
    except Exception:
        # 其他异常（如崩溃），标记为失败并回退
        _mark_amd_smi_failed()
    return gpu_map


def _map_amd_gpu_via_rocm_smi() -> dict:
    """使用 rocm-smi 获取 GPU 设备信息 (备选方案)"""
    gpu_map = {}
    try:
        # 获取产品名称
        product_json = _run_rocm_smi(["--showproductname", "--json"])
        # 获取显存信息 (VRAM)
        mem_json = _run_rocm_smi(["--showmeminfo", "vram", "--json"])

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


def _run_amd_smi_metric() -> List[dict]:
    """运行 amd-smi metric --json 并返回解析后的结果"""
    global _AMD_SMI_WORKING
    if not _AMD_SMI_WORKING:
        return []
    try:
        # 添加超时防止命令卡死，timeout 设为 5 秒
        result = subprocess.run(
            ["amd-smi", "metric", "--json"],
            capture_output=True,
            check=True,
            text=True,
            timeout=5,
        )
        return json.loads(result.stdout)
    except subprocess.TimeoutExpired:
        # 命令超时，标记为失败并回退
        _mark_amd_smi_failed()
        return []
    except (subprocess.CalledProcessError, json.JSONDecodeError):
        # 命令执行失败，标记为失败并回退
        _mark_amd_smi_failed()
        return []
    except Exception:
        # 其他异常（如崩溃），标记为失败并回退
        _mark_amd_smi_failed()
        return []


def _run_rocm_smi(args: list) -> dict:
    """
    辅助函数：运行 rocm-smi 并返回解析后的 JSON
    """
    try:
        cmd = ["rocm-smi"] + args
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return json.loads(result.stdout)
    except (subprocess.CalledProcessError, json.JSONDecodeError):
        return {}


class AMDCollector(H):
    def __init__(self, gpu_map, max_mem_value, system: str = "Linux"):
        super().__init__()
        self.gpu_map = gpu_map
        self.max_mem_value = max_mem_value
        self._system = system
        self._use_amd_smi = system == "Linux" and _USE_AMD_SMI

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

        # GPU Power (W) - Windows 下 hipinfo 不支持功耗监控
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
            # Linux 使用 amd-smi/rocm-smi
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

    def _parse_metric_value(self, value) -> Optional[float]:
        """解析指标值，处理 'N/A' 字符串等情况"""
        if value is None or value == "N/A":
            return None
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, dict):
            # 处理带单位的值，如 {"value": 29, "unit": "C"}
            val = value.get("value")
            if val is not None and val != "N/A":
                try:
                    return float(val)
                except (ValueError, TypeError):
                    pass
        if isinstance(value, str):
            try:
                return float(value)
            except (ValueError, TypeError):
                pass
        return None

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

    def get_utilization_usage(self) -> dict:
        """
        获取指定 GPU 设备的利用率
        amd-smi: usage 字段 (可能为 "N/A")
        rocm-smi: "GPU use (%)"
        Windows: hipinfo 不支持 GPU 利用率监控
        """
        util_infos = {}

        if self._use_amd_smi:
            metrics = _run_amd_smi_metric()
            for device in metrics:
                gpu_id = str(device.get("gpu", ""))
                util_infos[gpu_id] = {
                    "key": self.util_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Utilization (%)",
                    "value": math.nan,
                    "config": self.per_util_configs.get(f"GPU {gpu_id}"),
                }
                # amd-smi 的 usage 字段通常为 "N/A"，尝试解析
                usage = device.get("usage")
                if usage is not None and usage != "N/A":
                    try:
                        util_infos[gpu_id]["value"] = float(usage)
                    except (ValueError, TypeError):
                        pass
        else:
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
        amd-smi: mem_usage.total_vram 和 mem_usage.used_vram 计算
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

        # Linux 系统
        if self._use_amd_smi:
            metrics = _run_amd_smi_metric()
            for device in metrics:
                gpu_id = str(device.get("gpu", ""))
                mem_infos[gpu_id] = {
                    "key": self.memory_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Memory Allocated (%)",
                    "value": math.nan,
                    "config": self.per_memory_configs.get(f"GPU {gpu_id}"),
                }

                mem_usage = device.get("mem_usage", {})
                total_vram = self._parse_metric_value(mem_usage.get("total_vram"))
                used_vram = self._parse_metric_value(mem_usage.get("used_vram"))

                if total_vram is not None and used_vram is not None and total_vram > 0:
                    mem_infos[gpu_id]["value"] = (used_vram / total_vram) * 100
        else:
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
        amd-smi: mem_usage.used_vram (单位已是 MB)
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

        # Linux 系统
        if self._use_amd_smi:
            metrics = _run_amd_smi_metric()
            for device in metrics:
                gpu_id = str(device.get("gpu", ""))
                mem_value_infos[gpu_id] = {
                    "key": self.mem_value_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Memory Allocated (MB)",
                    "value": math.nan,
                    "config": self.per_mem_value_configs.get(f"GPU {gpu_id}"),
                }

                mem_usage = device.get("mem_usage", {})
                used_vram = self._parse_metric_value(mem_usage.get("used_vram"))
                if used_vram is not None:
                    mem_value_infos[gpu_id]["value"] = used_vram
        else:
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
        amd-smi: temperature.edge.value
        rocm-smi: "Temperature (Sensor junction) (C)" 或 "Temperature (Sensor edge) (C)"
        Windows: hipinfo 不支持温度监控
        """
        temp_infos = {}

        if self._use_amd_smi:
            metrics = _run_amd_smi_metric()
            for device in metrics:
                gpu_id = str(device.get("gpu", ""))

                temp_infos[gpu_id] = {
                    "key": self.temp_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Temperature (°C)",
                    "value": math.nan,
                    "config": self.per_temp_configs.get(f"GPU {gpu_id}"),
                }

                # 优先取 hotspot 温度，其次取 edge 温度
                temp_data = device.get("temperature", {})
                temp = self._parse_metric_value(temp_data.get("hotspot")) or self._parse_metric_value(temp_data.get("edge"))
                if temp is not None:
                    temp_infos[gpu_id]["value"] = temp
        else:
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
        amd-smi: power.socket_power (可能为 "N/A")
        rocm-smi: "Average Graphics Package Power (W)" 或 "Socket Power (W)"
        Windows: hipinfo 不支持功耗监控
        """
        power_infos = {}

        if self._use_amd_smi:
            metrics = _run_amd_smi_metric()
            for device in metrics:
                gpu_id = str(device.get("gpu", ""))

                power_infos[gpu_id] = {
                    "key": self.power_key.format(gpu_index=gpu_id),
                    "name": f"GPU {gpu_id} Power (W)",
                    "value": math.nan,
                    "config": self.per_power_configs.get(f"GPU {gpu_id}"),
                }

                power_data = device.get("power", {})
                socket_power = self._parse_metric_value(power_data.get("socket_power"))
                if socket_power is not None:
                    power_infos[gpu_id]["value"] = socket_power
        else:
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