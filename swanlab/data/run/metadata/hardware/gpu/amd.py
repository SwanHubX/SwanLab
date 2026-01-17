"""
@author: KMnO4-zx
@file: amd.py
@time: 2025-12-25
@description: AMD ROCm GPU 信息采集 (稳定性优化版)
"""

import json
import math
import platform
import subprocess
import re
import os
import glob
from typing import Optional, Tuple

from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index


def get_amd_gpu_info() -> HardwareFuncResult:
    """
    获取 AMD GPU 信息
    """
    system_name = platform.system()
    
    if system_name not in ["Linux", "Windows"]:
        return None, None

    info = {"driver": None, "gpu": None}
    collector = None
    
    try:
        if system_name == "Linux":
            driver, gpu_map = map_amd_gpu_linux()
        else:
            driver, gpu_map = map_amd_gpu_windows()
            
        info["driver"] = driver
        info["gpu"] = gpu_map
        
        # 计算最大显存用于图表配置
        max_mem_value = 0
        for gpu_id in gpu_map:
            try:
                mem_str = gpu_map[gpu_id]["memory"]
                digits = re.findall(r"([\d\.]+)", mem_str)
                if digits:
                    mem_val = float(digits[0])
                    max_mem_value = max(max_mem_value, mem_val)
            except (ValueError, TypeError):
                continue
        
        max_mem_value = int(max_mem_value * 1024)
        
        collector = AMDCollector(gpu_map, max_mem_value, system_name)
        
    except Exception:
        if all(v is None for v in info.values()):
            return None, None
            
    return info, collector


# ==========================================
# Linux Logic (Initialization)
# ==========================================

def map_amd_gpu_linux() -> Tuple[Optional[str], dict]:
    """
    初始化阶段：尝试获取显卡列表。
    虽然 rocm-smi 可能不稳定，但我们只在启动时跑一次。
    如果 rocm-smi 失败，尝试通过 lspci 或 sysfs 兜底。
    """
    driver_version = "Unknown"
    
    # 1. 尝试获取版本
    try:
        # 读取 /sys/module/amdgpu/version (最安全的方式)
        if os.path.exists("/sys/module/amdgpu/version"):
            with open("/sys/module/amdgpu/version", "r") as f:
                driver_version = f.read().strip()
        else:
            # 回退到命令
            res = subprocess.run(["rocm-smi", "--showdriverversion", "--json"], capture_output=True, text=True)
            if "Driver version" in res.stdout:
                driver_version = json.loads(res.stdout).get("card0", {}).get("Driver version", "Unknown")
    except Exception:
        pass

    gpu_map = {}
    
    # 2. 尝试获取设备列表
    # 优先使用 sysfs 扫描，避免调用 rocm-smi 导致崩溃
    try:
        cards = glob.glob("/sys/class/drm/card*")
        # 过滤出是 GPU 的卡（排除 card0 可能是集成显卡的情况，如果都是 AMD，通常是 card0, card1...）
        # 这里简单假设所有 /sys/class/drm/cardX 且包含 device/vendor 的都是目标卡
        
        amd_cards = []
        for card_path in cards:
            # 排除 renderD* 或 controlD*
            if not re.search(r"card\d+$", card_path):
                continue
                
            # 检查 vendor id, AMD 是 0x1002
            vendor_path = os.path.join(card_path, "device/vendor")
            if os.path.exists(vendor_path):
                with open(vendor_path, "r") as f:
                    vendor = f.read().strip()
                if "0x1002" in vendor:
                    idx = card_path.split("card")[-1]
                    amd_cards.append(idx)
        
        amd_cards.sort(key=lambda x: int(x))
        
        for idx in amd_cards:
            # 显存大小读取
            mem_total_path = f"/sys/class/drm/card{idx}/device/mem_info_vram_total"
            mem_str = "0GB"
            if os.path.exists(mem_total_path):
                with open(mem_total_path, "r") as f:
                    mem_bytes = int(f.read().strip())
                    mem_gb = round(mem_bytes / (1024**3))
                    mem_str = f"{mem_gb}GB"
            
            # 名字读取 (sysfs 中没有漂亮的 Marketing Name，只有 device id)
            # 为了获取名字，我们可以尝试 rocm-smi 一次，如果失败就用 "AMD GPU"
            name = "AMD GPU" 
            
            gpu_map[str(idx)] = {"name": name, "memory": mem_str}

    except Exception:
        # 如果 Sysfs 失败，最后尝试一次 rocm-smi
        try:
            data = json.loads(subprocess.run(["rocm-smi", "--showproductname", "--json"], capture_output=True, text=True).stdout)
            for k, v in data.items():
                idx = k.replace("card", "")
                gpu_map[idx] = {"name": v.get("Card Series", "AMD GPU"), "memory": "0GB"} # 显存稍后补
        except Exception:
            pass

    if len(gpu_map) == 0:
        return (None, None)

    return driver_version, gpu_map


# ==========================================
# Windows Logic (hipinfo)
# ==========================================

def map_amd_gpu_windows() -> Tuple[Optional[str], dict]:
    driver_version = "Unknown" 
    gpu_map = {}
    try:
        output = subprocess.run(["hipinfo"], capture_output=True, text=True, check=True).stdout
        devices = re.split(r"device#\s+(\d+)", output)
        if len(devices) < 2:
            gpu_map = _parse_hipinfo_text_block_windows("0", output)
        else:
            for i in range(1, len(devices), 2):
                dev_id = devices[i].strip()
                block_content = devices[i+1]
                gpu_map.update(_parse_hipinfo_text_block_windows(dev_id, block_content))
    except Exception:
        return None, None
    return driver_version, gpu_map

def _parse_hipinfo_text_block_windows(dev_id: str, content: str) -> dict:
    name = "AMD GPU"
    memory = "0GB"
    name_match = re.search(r"Name:\s*(.+)", content)
    if name_match: name = name_match.group(1).strip()
    mem_match = re.search(r"totalGlobalMem:\s*([\d\.]+)\s*GB", content)
    if mem_match: memory = f"{mem_match.group(1)}GB"
    return {dev_id: {"name": name, "memory": memory}}


# ==========================================
# Collector Class
# ==========================================

class AMDCollector(H):
    def __init__(self, gpu_map, max_mem_value, system_name):
        super().__init__()
        self.gpu_map = gpu_map
        self.max_mem_value = max_mem_value
        self.system_name = system_name

        self.util_key = generate_key("gpu.{gpu_index}.pct")
        self.memory_key = generate_key("gpu.{gpu_index}.mem.pct")
        self.mem_value_key = generate_key("gpu.{gpu_index}.mem.value")
        self.temp_key = generate_key("gpu.{gpu_index}.temp")
        self.power_key = generate_key("gpu.{gpu_index}.power")

        util_config = HardwareConfig(y_range=(0, 100), chart_index=random_index(), chart_name="GPU Utilization (%)")
        mem_pct_config = HardwareConfig(y_range=(0, 100), chart_index=random_index(), chart_name="GPU Memory Allocated (%)")
        mem_val_config = HardwareConfig(y_range=(0, self.max_mem_value), chart_index=random_index(), chart_name="GPU Memory Allocated (MB)")
        temp_config = HardwareConfig(y_range=(0, None), chart_index=random_index(), chart_name="GPU Temperature (°C)")
        power_config = HardwareConfig(y_range=(0, None), chart_index=random_index(), chart_name="GPU Power (W)")

        self.per_util_configs = {}
        self.per_memory_configs = {}
        self.per_mem_value_configs = {}
        self.per_temp_configs = {}
        self.per_power_configs = {}

        for gpu_id in self.gpu_map:
            metric_name = f"GPU {gpu_id}"
            self.per_util_configs[metric_name] = util_config.clone(metric_name=metric_name)
            self.per_memory_configs[metric_name] = mem_pct_config.clone(metric_name=metric_name)
            self.per_mem_value_configs[metric_name] = mem_val_config.clone(metric_name=metric_name)
            self.per_temp_configs[metric_name] = temp_config.clone(metric_name=metric_name)
            self.per_power_configs[metric_name] = power_config.clone(metric_name=metric_name)

    def collect(self) -> HardwareInfoList:
        if self.system_name == "Linux":
            return self._collect_linux_sysfs()
        else:
            return self._collect_windows()

    def _collect_linux_sysfs(self) -> HardwareInfoList:
        """
        Linux 采集 (Sysfs 版) - 不执行任何命令，只读文件，解决 Core Dump。
        """
        result = []
        
        for gpu_id in self.gpu_map:
            # 路径构造: /sys/class/drm/card0/device/
            base_path = f"/sys/class/drm/card{gpu_id}/device"
            
            # 1. Utilization
            # 文件通常是 gpu_busy_percent，内容是 0-100 的整数
            util_val = math.nan
            try:
                with open(os.path.join(base_path, "gpu_busy_percent"), "r") as f:
                    util_val = float(f.read().strip())
            except Exception:
                pass
            result.append(self._build_metric(gpu_id, self.util_key, self.per_util_configs, util_val))

            # 2. Memory
            # mem_info_vram_used (Bytes) / mem_info_vram_total (Bytes)
            mem_pct = math.nan
            mem_used_mb = math.nan
            try:
                with open(os.path.join(base_path, "mem_info_vram_used"), "r") as f:
                    used_bytes = int(f.read().strip())
                with open(os.path.join(base_path, "mem_info_vram_total"), "r") as f:
                    total_bytes = int(f.read().strip())
                
                if total_bytes > 0:
                    mem_pct = (used_bytes / total_bytes) * 100
                    mem_used_mb = used_bytes / (1024 * 1024)
            except Exception:
                pass
            
            result.append(self._build_metric(gpu_id, self.memory_key, self.per_memory_configs, mem_pct))
            result.append(self._build_metric(gpu_id, self.mem_value_key, self.per_mem_value_configs, mem_used_mb))

            # 3. Temp & Power
            # 位于 hwmon 子目录，例如 /sys/class/drm/card0/device/hwmon/hwmon1/
            # 需要动态寻找 hwmon 目录
            temp_val = math.nan
            power_val = math.nan
            
            try:
                hwmon_pattern = os.path.join(base_path, "hwmon", "hwmon*")
                hwmons = glob.glob(hwmon_pattern)
                if hwmons:
                    hwmon_path = hwmons[0] # 通常只有一个
                    
                    # Temp: temp1_input (millidegrees C) -> / 1000
                    # 某些卡可能是 edge_input 或 junction_input，通常 temp1_input 是主要温度
                    if os.path.exists(os.path.join(hwmon_path, "temp1_input")):
                        with open(os.path.join(hwmon_path, "temp1_input"), "r") as f:
                            temp_val = float(f.read().strip()) / 1000.0
                            
                    # Power: power1_average (microwatts) -> / 1000000 -> W
                    # 有些是 power1_input
                    p_file = os.path.join(hwmon_path, "power1_average")
                    if not os.path.exists(p_file):
                        p_file = os.path.join(hwmon_path, "power1_input")
                        
                    if os.path.exists(p_file):
                        with open(p_file, "r") as f:
                            power_val = float(f.read().strip()) / 1000000.0
            except Exception:
                pass

            result.append(self._build_metric(gpu_id, self.temp_key, self.per_temp_configs, temp_val))
            result.append(self._build_metric(gpu_id, self.power_key, self.per_power_configs, power_val))

        return result

    def _collect_windows(self) -> HardwareInfoList:
        result = []
        output = ""
        try:
            output = subprocess.run(["hipinfo"], capture_output=True, text=True).stdout
        except Exception:
            pass

        devices_content = {}
        parts = re.split(r"device#\s+(\d+)", output)
        if len(parts) >= 2:
            for i in range(1, len(parts), 2):
                devices_content[parts[i].strip()] = parts[i+1]
        else:
            devices_content["0"] = output

        for gpu_id in self.gpu_map:
            content = devices_content.get(gpu_id, "")
            
            mem_total = self._extract_value_windows(content, r"memInfo\.total:\s*([\d\.]+)\s*GB")
            mem_free = self._extract_value_windows(content, r"memInfo\.free:\s*([\d\.]+)\s*GB")
            
            mem_pct = math.nan
            mem_used_mb = math.nan
            
            if mem_total is not None and mem_free is not None and mem_total > 0:
                mem_used = mem_total - mem_free
                mem_used_mb = mem_used * 1024
                mem_pct = (mem_used / mem_total) * 100
                
            result.append(self._build_metric(gpu_id, self.memory_key, self.per_memory_configs, mem_pct))
            result.append(self._build_metric(gpu_id, self.mem_value_key, self.per_mem_value_configs, mem_used_mb))
            result.append(self._build_metric(gpu_id, self.temp_key, self.per_temp_configs, math.nan))
            result.append(self._build_metric(gpu_id, self.util_key, self.per_util_configs, math.nan))
            result.append(self._build_metric(gpu_id, self.power_key, self.per_power_configs, math.nan))
            
        return result

    def _extract_value_windows(self, content: str, regex: str) -> Optional[float]:
        try:
            match = re.search(regex, content, re.IGNORECASE)
            if match: return float(match.group(1))
        except Exception: pass
        return None

    def _build_metric(self, gpu_id, key_template, config_map, value):
        key = key_template.format(gpu_index=gpu_id)
        final_val = math.nan
        if value is not None and not math.isnan(value):
            try: final_val = float(value)
            except Exception: pass
        return {
            "key": key,
            "name": config_map[f"GPU {gpu_id}"].chart_name,
            "value": final_val,
            "config": config_map[f"GPU {gpu_id}"],
        }