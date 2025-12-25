"""
@author: ZhiningMa
@file: iluvatar.py
@time: 2025-12-24
@description: Iluvatar GPU 信息采集
"""
import platform
import subprocess
from typing import Optional, Tuple
from swanlab.log import swanlog
from ..type import HardwareCollector as H
from ..type import HardwareConfig, HardwareFuncResult, HardwareInfoList
from ..utils import generate_key, random_index

def get_iluvatar_gpu_info() ->HardwareFuncResult:
    """
    获取Iluvatar GPU信息，包括驱动版本、设备信息等
    """
    # Iluvatar GPU只支持Linux系统
    if platform.system() != "Linux":
        return None, None    
    
    info = {"driver": None,  "gpu": None}
    collector = None
    try:
        driver, gpu_map = map_iluvatar_gpu()
        info["driver"] = driver
        info["gpu"] = gpu_map
        max_mem_value = 0
        for gpu_id in gpu_map:
            mem_value = int(gpu_map[gpu_id]["memory"])
            max_mem_value = max(max_mem_value, mem_value)
        max_mem_value *= 1024
        collector = IluvatarCollector(gpu_map, max_mem_value)
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None, None
    return info, collector
    
    
def map_iluvatar_gpu() -> Tuple[Optional[str], dict]:
    """
    获取 Iluvatar GPU信息，包括驱动版本、设备信息等，例如：
    driver: '4.4.0'
    gpu_map: {"0": { "name": "Iluvatar BI-V150", "memory": "32"}, "1": { "name": "Iluvatar BI-V150", "memory": "32"}, ...}
    """

    # 运行 ixsmi 命令获取 GPU 信息，添加 index 字段
    output_str = subprocess.run(
        ["ixsmi", "--query-gpu=index,driver_version,name,memory.total", "--format=csv"],
        capture_output=True, 
        check=True, 
        text=True,
        encoding='utf-8'
    ).stdout
    
    # 分割输出行
    lines = output_str.strip().split('\n')
    
    if len(lines) < 2:
        return None, {}
    
    driver = None
    gpu_map = {}
    
    # 处理每一行数据（跳过标题行）
    for line in lines[1:]:
        # 分割 CSV 行，处理可能包含空格的情况
        parts = [part.strip() for part in line.split(',')]
        
        if len(parts) >= 4:
            # 获取 GPU ID（使用 index 字段）
            gpu_id = parts[0]
            
            # 驱动版本（只取第一个GPU的驱动版本）
            if driver is None:
                driver = parts[1]
            
            # GPU 名称
            gpu_name = parts[2]
            
            # 内存处理
            gpu_memory_str = parts[3]
            if "MiB" in gpu_memory_str:
                # 提取数字并转换为GB
                memory_mib = int(gpu_memory_str.replace(" MiB", "").strip())
                memory_gb = str(memory_mib // 1024)
            elif "GiB" in gpu_memory_str:
                # 如果是GiB，直接提取数字
                memory_gb = str(int(gpu_memory_str.replace(" GiB", "").strip()))
            else:
                # 如果不是标准格式，保持原样
                memory_gb = gpu_memory_str
            
            gpu_map[gpu_id] = {
                "name": gpu_name,
                "memory": memory_gb,
            }
    
    return driver, gpu_map
        


class IluvatarCollector(H):
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
        self.memory_value_key = generate_key("gpu.{gpu_index}.mem.value")
        memory_value_config = HardwareConfig(
            y_range=(0, self.max_mem_value),
            chart_index=random_index(),
            chart_name="GPU Memory Allocated (MB)",
        )
        self.per_memory_value_configs = {}

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
            self.per_memory_value_configs[metric_name] = memory_value_config.clone(metric_name=metric_name)
            self.per_temp_configs[metric_name] = temp_config.clone(metric_name=metric_name)
            self.per_power_configs[metric_name] = power_config.clone(metric_name=metric_name)
        
    def collect(self) -> HardwareInfoList:
        result: HardwareInfoList = []
        usage_methods = [
            self.get_utilization_usage,
            self.get_memory_rate_usage,
            self.get_memory_value_usage,
            self.get_temperature_usage,
            self.get_power_usage,
        ]

        for method in usage_methods:
            result.extend(method().values())
        return result
     
    def _run_ixsmi_query(self, query_type: str) -> list:
        """
        运行 ixsmi 查询并解析结果
        query_type: utilization.gpu, utilization.memory, temperature.gpu, memory.used, board.power.draw
        """

        output_str = subprocess.run(
            ["ixsmi", f"--query-gpu={query_type}", "--format=csv"],
            capture_output=True,
            check=True,
            text=True,
            encoding='utf-8'
        ).stdout
        
        # 分割输出行，跳过标题行
        lines = output_str.strip().split('\n')
        if len(lines) < 2:
            return []
        
        # 解析每一行数据
        results = []
        for line in lines[1:]:  # 跳过标题行
            parts = line.strip().split(',')
            if len(parts) > 0:
                # 清理数据：移除单位
                value_str = parts[0].strip()
                
                # 根据不同的查询类型处理数据
                if query_type == "utilization.gpu":
                    # 移除 % 符号
                    value = float(value_str.replace('%', '').strip())
                elif query_type == "utilization.memory":
                    # 移除 % 符号
                    value = float(value_str.replace('%', '').strip())
                elif query_type == "temperature.gpu":
                    # 移除 C 符号
                    value = float(value_str.replace('C', '').strip())
                elif query_type == "memory.used":
                    # 移除 MiB 符号并转换为MB
                    value = float(value_str.replace('MiB', '').strip())
                elif query_type == "board.power.draw":
                    # 移除 W 符号
                    value = float(value_str.replace('W', '').strip())
                else:
                    swanlog.warning(f"Unknown query type in _run_ixsmi_query: {query_type}")
                    value = 0.0
                
                results.append(value)
        
        return results
            

    
    def get_utilization_usage(self) -> dict:
        """
        获取所有GPU设备的利用率（%）
        """
        gpu_utilizations = self._run_ixsmi_query("utilization.gpu")
        usage_infos = {}
        
        for i, gpu_usage in enumerate(gpu_utilizations):
            usage_infos[i] = {
                "key": self.util_key.format(gpu_index=i),
                "name": f"GPU {i} Utilization (%)",
                "value": gpu_usage,
                "config": self.per_util_configs[f"GPU {i}"],
            }
        
        return usage_infos
    
    def get_memory_rate_usage(self) -> dict:
        """
        获取所有GPU设备的内存占用率（%）
        """
        memory_utilizations = self._run_ixsmi_query("utilization.memory")
        mem_infos = {}
        
        for i, mem_usage_rate in enumerate(memory_utilizations):
            mem_infos[i] = {
                "key": self.memory_key.format(gpu_index=i),
                "name": f"GPU {i} Memory Allocated (%)",
                "value": mem_usage_rate,
                "config": self.per_memory_configs[f"GPU {i}"],
            }
        
        return mem_infos
    
    def get_memory_value_usage(self) -> dict:
        """
        获取所有GPU设备的内存占用值（MB）
        """
        memory_used = self._run_ixsmi_query("memory.used")
        mem_value_infos = {}
        
        for i, mem_usage_value in enumerate(memory_used):
            mem_value_infos[i] = {
                "key": self.memory_value_key.format(gpu_index=i),
                "name": f"GPU {i} Memory Allocated (MB)",
                "value": mem_usage_value,
                "config": self.per_memory_value_configs[f"GPU {i}"],
            }
        
        return mem_value_infos
    
    def get_temperature_usage(self) -> dict:
        """
        获取所有GPU设备的温度（°C）
        """
        temperatures = self._run_ixsmi_query("temperature.gpu")
        temp_infos = {}
        
        for i, temp in enumerate(temperatures):
            temp_infos[i] = {
                "key": self.temp_key.format(gpu_index=i),
                "name": f"GPU {i} Temperature (°C)",
                "value": temp,
                "config": self.per_temp_configs[f"GPU {i}"],
            }
        
        return temp_infos
    
    def get_power_usage(self) -> dict:
        """
        获取所有GPU设备的功耗（W）
        """
        power_draws = self._run_ixsmi_query("board.power.draw")
        power_infos = {}
        
        for i, power in enumerate(power_draws):
            power_infos[i] = {
                "key": self.power_key.format(gpu_index=i),
                "name": f"GPU {i} Power (W)",
                "value": power,
                "config": self.per_power_configs[f"GPU {i}"],
            }
        
        return power_infos
