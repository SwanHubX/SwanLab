"""
@author: KMnO4-zx
@file: amd.py
@time: 2025-12-24
@description: AMD ROCm GPU 信息采集 (稳定性优化版)
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
    获取 AMD GPU 信息，自动判断 Linux (rocm-smi) 或 Windows (hipinfo)
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
                # 提取数字 (支持 107.87 GB 这种格式)
                digits = re.findall(r"([\d\.]+)", mem_str)
                if digits:
                    mem_val = float(digits[0])
                    max_mem_value = max(max_mem_value, mem_val)
            except (ValueError, TypeError):
                continue
        
        # 将 GB 转为 MB (保留精度后转整型，防止 16.5GB 被截断)
        max_mem_value = int(max_mem_value * 1024)
        
        collector = AMDCollector(gpu_map, max_mem_value, system_name)
        
    except Exception:
        if all(v is None for v in info.values()):
            return None, None
            
    return info, collector


# ==========================================
# Linux Logic (rocm-smi)
# ==========================================

def map_amd_gpu_linux() -> Tuple[Optional[str], dict]:
    driver_version = "Unknown"
    try:
        version_out = subprocess.run(["rocm-smi", "--showdriverversion", "--json"], capture_output=True, text=True).stdout
        # 如果 rocm-smi 获取不到，尝试 amd-smi
        if "Driver version" not in version_out:
             res = subprocess.run(["amd-smi", "version"], capture_output=True, text=True)
             match = re.search(r"ROCm version:\s*([\d\.]+)", res.stdout)
             if match: driver_version = match.group(1)
    except Exception:
        pass

    gpu_map = {}
    try:
        # 合并查询
        data = _run_rocm_smi_linux(["--showproductname", "--showmeminfo", "vram", "--json"])
        for card_key, info in data.items():
            idx = card_key.replace("card", "")
            name = info.get("Card Series", "AMD GPU")
            mem_bytes_str = info.get("VRAM Total Memory (B)", "0")
            try:
                mem_gb = round(int(mem_bytes_str) / (1024**3))
                mem_str = f"{mem_gb}GB"
            except:
                mem_str = "0GB"
            gpu_map[str(idx)] = {"name": name, "memory": mem_str}
    except Exception:
        return driver_version, {}

    return driver_version, gpu_map

def _run_rocm_smi_linux(args: list) -> dict:
    try:
        cmd = ["rocm-smi"] + args
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=5)
        return json.loads(result.stdout)
    except Exception:
        return {}


# ==========================================
# Windows Logic (hipinfo)
# ==========================================

def map_amd_gpu_windows() -> Tuple[Optional[str], dict]:
    """
    Windows 解析 hipinfo 输出
    """
    driver_version = "Unknown" # hipinfo 输出中似乎没有明确的 Driver Version，暂定 Unknown
    gpu_map = {}
    
    try:
        output = subprocess.run(["hipinfo"], capture_output=True, text=True, check=True).stdout
        
        # 按 device# 分割，第一部分通常为空或头部信息
        # 输出示例: device# 0 ...
        devices = re.split(r"device#\s+(\d+)", output)
        
        if len(devices) < 2:
            # 只有一个设备或格式不匹配，尝试直接全量解析
            gpu_map = _parse_hipinfo_text_block_windows("0", output)
        else:
            # devices[0] 是头部，devices[1] 是 id，devices[2] 是内容...
            for i in range(1, len(devices), 2):
                dev_id = devices[i].strip()
                block_content = devices[i+1]
                gpu_map.update(_parse_hipinfo_text_block_windows(dev_id, block_content))
                
    except Exception:
        return None, None
        
    return driver_version, gpu_map

def _parse_hipinfo_text_block_windows(dev_id: str, content: str) -> dict:
    """
    解析 Windows hipinfo 的单个设备块
    Name: AMD Radeon(TM) 8060S Graphics
    totalGlobalMem: 107.87 GB
    """
    name = "AMD GPU"
    memory = "0GB"
    
    # 解析名称
    name_match = re.search(r"Name:\s*(.+)", content)
    if name_match:
        name = name_match.group(1).strip()
        
    # 解析总显存
    mem_match = re.search(r"totalGlobalMem:\s*([\d\.]+)\s*GB", content)
    if mem_match:
        # 直接使用获取到的浮点数显存，保留原貌
        memory = f"{mem_match.group(1)}GB"
        
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
            return self._collect_linux()
        else:
            return self._collect_windows()

    def _collect_linux(self) -> HardwareInfoList:
        try:
            data = _run_rocm_smi_linux(["--showuse", "--showmemuse", "--showtemp", "--showpower", "--json"])
        except Exception:
            return []

        result = []
        for gpu_id in self.gpu_map:
            card_key = f"card{gpu_id}"
            info = data.get(card_key, {})
            
            # Util
            u_val = info.get("GPU use (%)")
            result.append(self._build_metric(gpu_id, self.util_key, self.per_util_configs, u_val))

            # Mem
            m_pct_val = info.get("GPU Memory use (%)") or info.get("VRAM use (%)")
            m_mb_val = None
            if m_pct_val is not None:
                try:
                    mem_str = self.gpu_map[gpu_id]["memory"]
                    total_mb = float(re.findall(r"[\d\.]+", mem_str)[0]) * 1024
                    m_mb_val = float(m_pct_val) * 0.01 * total_mb
                except: pass
            result.append(self._build_metric(gpu_id, self.memory_key, self.per_memory_configs, m_pct_val))
            result.append(self._build_metric(gpu_id, self.mem_value_key, self.per_mem_value_configs, m_mb_val))

            # Temp
            t_val = info.get("Temperature (Sensor junction) (C)") or info.get("Temperature (Sensor edge) (C)")
            result.append(self._build_metric(gpu_id, self.temp_key, self.per_temp_configs, t_val))

            # Power
            p_val = info.get("Average Graphics Package Power (W)") or info.get("Socket Power (W)")
            if p_val is None:
                for k, v in info.items():
                    if "Power" in k and "(W)" in k: p_val = v; break
            result.append(self._build_metric(gpu_id, self.power_key, self.per_power_configs, p_val))
        return result

    def _collect_windows(self) -> HardwareInfoList:
        """
        Windows 采集逻辑：
        仅计算显存 (Total - Free)。
        Temp, Power, Util 设为 NaN。
        """
        result = []
        output = ""
        try:
            output = subprocess.run(["hipinfo"], capture_output=True, text=True).stdout
        except:
            pass

        # 如果有多设备，这里其实需要按 device# 分割，但为简化起见，我们假设 output 里
        # 包含了所有信息，且我们使用全文本正则搜索（不够严谨但能跑通单卡）
        # 严谨做法是再次调用 map 里的分割逻辑
        
        devices_content = {}
        parts = re.split(r"device#\s+(\d+)", output)
        if len(parts) >= 2:
            for i in range(1, len(parts), 2):
                devices_content[parts[i].strip()] = parts[i+1]
        else:
            devices_content["0"] = output

        for gpu_id in self.gpu_map:
            content = devices_content.get(gpu_id, "")
            
            # 1. Memory Calculation (这是 Windows 唯一能动的)
            # memInfo.total: 107.87 GB
            # memInfo.free: 107.72 GB (100%)
            mem_total = self._extract_value_windows(content, r"memInfo\.total:\s*([\d\.]+)\s*GB")
            mem_free = self._extract_value_windows(content, r"memInfo\.free:\s*([\d\.]+)\s*GB")
            
            mem_pct = math.nan
            mem_used_mb = math.nan
            
            if mem_total is not None and mem_free is not None and mem_total > 0:
                mem_used = mem_total - mem_free
                mem_used_mb = mem_used * 1024 # GB to MB
                mem_pct = (mem_used / mem_total) * 100
                
            result.append(self._build_metric(gpu_id, self.memory_key, self.per_memory_configs, mem_pct))
            result.append(self._build_metric(gpu_id, self.mem_value_key, self.per_mem_value_configs, mem_used_mb))
            
            # 2. Others (Temp, Util, Power) -> NaN
            result.append(self._build_metric(gpu_id, self.temp_key, self.per_temp_configs, math.nan))
            result.append(self._build_metric(gpu_id, self.util_key, self.per_util_configs, math.nan))
            result.append(self._build_metric(gpu_id, self.power_key, self.per_power_configs, math.nan))
            
        return result

    def _extract_value_windows(self, content: str, regex: str) -> Optional[float]:
        try:
            match = re.search(regex, content, re.IGNORECASE)
            if match: return float(match.group(1))
        except: pass
        return None

    def _build_metric(self, gpu_id, key_template, config_map, value):
        key = key_template.format(gpu_index=gpu_id)
        final_val = math.nan
        if value is not None:
            try: final_val = float(value)
            except: pass
        return {
            "key": key,
            "name": config_map[f"GPU {gpu_id}"].chart_name,
            "value": final_val,
            "config": config_map[f"GPU {gpu_id}"],
        }