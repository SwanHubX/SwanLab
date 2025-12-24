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
from typing import Optional, Tuple

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


def get_amd_gpu_info() -> HardwareFuncResult:
    """
    获取 AMD GPU 信息，包括驱动版本、设备信息等
    """
    # AMD ROCm 主要运行在 Linux 系统
    if platform.system() != "Linux":
        return None, None

    info = {"driver": None, "gpu": None}
    collector = None
    try:
        driver, gpu_map = map_amd_gpu()
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
        collector = AMDCollector(gpu_map, max_mem_value)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector


def map_amd_gpu() -> Tuple[Optional[str], dict]:
    """
    获取 AMD GPU 信息，包括 ROCm 版本、设备列表
    返回示例:
    ROCm Version: 6.4.3
    {'0': {'name': 'AMD Radeon Graphics', 'memory': '64GB'}, ...}
    """
    # 1. 获取版本信息 (优先使用 amd-smi version，因为它包含详细的 ROCm 版本)
    driver_version = None
    try:
        version_output = subprocess.run(
            ["amd-smi", "version"], capture_output=True, check=True, text=True
        ).stdout
        # 解析 "ROCm version: 6.4.3"
        match = re.search(r"ROCm version:\s*([\d\.]+)", version_output)
        if match:
            driver_version = match.group(1)
        else:
            # 备选：尝试从 rocm-smi 获取
            driver_version = "Unknown"
    except Exception:
        driver_version = "Unknown"

    # 2. 获取设备名称和显存信息 (使用 rocm-smi --json)
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
        return driver_version, {}

    return driver_version, gpu_map


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
    def __init__(self, gpu_map, max_mem_value):
        super().__init__()
        self.gpu_map = gpu_map
        self.max_mem_value = max_mem_value

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
        rocm-smi output key: "GPU use (%)"
        """
        output_json = _run_rocm_smi(["--showuse", "--json"])
        util_infos = {}
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")
            
            util_infos[gpu_id] = {
                "key": self.util_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Utilization (%)",
                "value": math.nan,
                "config": self.per_util_configs.get(f"GPU {gpu_id}"),
            }
            
            # 兼容不同版本的 Key
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
        rocm-smi output key: "GPU Memory use (%)" (不同版本可能略有差异，也可能是 "VRAM use (%)")
        """
        output_json = _run_rocm_smi(["--showmemuse", "--json"])
        mem_infos = {}

        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")
            mem_infos[gpu_id] = {
                "key": self.memory_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Memory Allocated (%)",
                "value": math.nan,
                "config": self.per_memory_configs.get(f"GPU {gpu_id}"),
            }

            # 尝试获取不同的可能的 Key
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
        由于 rocm-smi 有时直接给百分比，这里通过总显存 * 百分比计算，或者直接获取已用显存
        """
        # 为了准确，我们这里复用 showmemuse 的百分比 * 总显存 (从 self.max_mem_value 估算)
        # 或者使用 --showmeminfo vram 再次获取实时 usage (如果支持)
        # 简单起见，这里复用百分比计算逻辑，保持和 hygon.py 一致的逻辑流
        
        # 重新获取一次百分比
        output_json = _run_rocm_smi(["--showmemuse", "--json"])
        mem_value_infos = {}

        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")
            mem_value_infos[gpu_id] = {
                "key": self.mem_value_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Memory Allocated (MB)",
                "value": math.nan,
                "config": self.per_mem_value_configs.get(f"GPU {gpu_id}"),
            }

            val = info.get("GPU Memory use (%)") or info.get("VRAM use (%)")
            
            # 注意：self.max_mem_value 是所有卡中最大的显存，这里如果卡显存不一致可能有误差
            # 这种实现假设了多卡显存一致，或者仅作为近似值。
            # 更精确的做法是保存每张卡的 max_mem
            
            if val is not None:
                try:
                    # 获取当前卡的显存上限 (从 map 中找)
                    # 注意：map 中的 memory 是 "64GB" 字符串
                    current_card_mem_str = self.gpu_map.get(gpu_id, {}).get("memory", "0GB")
                    current_card_mem_mb = int(float(re.sub(r"[^0-9.]", "", current_card_mem_str))) * 1024
                    
                    mem_value_infos[gpu_id]["value"] = float(val) * 0.01 * current_card_mem_mb
                except (ValueError, TypeError):
                    pass
        return mem_value_infos

    def get_temperature_usage(self) -> dict:
        """
        获取指定 GPU 设备的温度(°C)
        rocm-smi output key 示例: "Temperature (Sensor edge) (C)"
        """
        output_json = _run_rocm_smi(["--showtemp", "--json"])
        temp_infos = {}
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")

            temp_infos[gpu_id] = {
                "key": self.temp_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Temperature (°C)",
                "value": math.nan,
                "config": self.per_temp_configs.get(f"GPU {gpu_id}"),
            }
            
            # 优先取 Junction 温度，其次取 Edge 温度 (和 Hygon 逻辑类似，Hygon 取了 Junction)
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
        rocm-smi output key 示例: "Average Graphics Package Power (W)" 或 "Socket Power (W)"
        """
        output_json = _run_rocm_smi(["--showpower", "--json"])
        power_infos = {}
        for card_key, info in output_json.items():
            gpu_id = card_key.replace("card", "")

            power_infos[gpu_id] = {
                "key": self.power_key.format(gpu_index=gpu_id),
                "name": f"GPU {gpu_id} Power (W)",
                "value": math.nan,
                "config": self.per_power_configs.get(f"GPU {gpu_id}"),
            }
            
            # 尝试不同的功耗 Key
            val = info.get("Average Graphics Package Power (W)") or info.get("Socket Power (W)")
            # 某些版本直接是 "Power (Socket)" (不带单位后缀但带括号) ? 根据用户 log 是 "Power (Socket)" 在表中
            # JSON 模式下 rocm-smi 通常会带全名。如果都不行，尝试模糊匹配
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