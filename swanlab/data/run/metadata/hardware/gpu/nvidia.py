"""
@author: cunyue
@file: nvidia.py
@time: 2024/12/3 20:12
@description: NVIDIA GPU信息采集
"""

import subprocess

import pynvml

from swanlab.log import swanlog
from ..type import HardwareFuncResult, HardwareCollector, HardwareInfoList, HardwareInfo
from ..utils import generate_key, HardwareConfig, random_index


def get_nvidia_gpu_info() -> HardwareFuncResult:
    """获取 GPU 信息"""

    info = {"driver": None, "cores": None, "type": [], "memory": [], "cuda": None}
    try:
        pynvml.nvmlInit()
    except Exception:  # noqa
        return None, None
    try:
        # 获取 NVIDIA 驱动版本信息
        nv_driver = pynvml.nvmlSystemGetDriverVersion()
        if isinstance(nv_driver, bytes):
            nv_driver = nv_driver.decode("utf-8")
        info["driver"] = nv_driver

        # 获取 CUDA 版本
        info["cuda"] = get_cuda_version()

        # 获取 NVIDIA GPU 数量
        info["cores"] = pynvml.nvmlDeviceGetCount()
        # 遍历每个 GPU，获取 GPU 信息
        for i in range(info["cores"]):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            # 获取 GPU 型号
            gpu_name = pynvml.nvmlDeviceGetName(handle)  # types: bytes | str
            if isinstance(gpu_name, bytes):  # Fix for pynvml 早期版本，关联 issue: #605
                gpu_name = gpu_name.decode("utf-8")
            info["type"].append(gpu_name)
            # 获取 GPU 的总显存, 单位为GB
            info["memory"].append(round(pynvml.nvmlDeviceGetMemoryInfo(handle).total / (1024**3)))
    except UnicodeDecodeError:  # 部分GPU型号无法解码
        return None, None
    except pynvml.NVMLError:
        pass
    finally:
        pynvml.nvmlShutdown()
    count = info["cores"]
    return info, None if not count else GpuCollector(count)


def get_cuda_version():
    """获取 CUDA 版本"""
    try:
        output = subprocess.check_output(["nvcc", "--version"]).decode("utf-8")
        for line in output.split('\n'):
            if "release" in line:
                version = line.split("release")[-1].strip().split(" ")[0][:-1]
                return version
    except Exception:  # noqa
        return None


class GpuCollector(HardwareCollector):

    def __init__(self, count: int):
        super().__init__()
        # GPU 利用率
        self.gpu_util_key = generate_key("gpu.{idx}.ptc")
        util_config = HardwareConfig(y_range=(0, 100), chart_name="GPU Utilization (%)", chart_index=random_index())
        # GPU 内存使用率
        self.gpu_mem_pct_key = generate_key("gpu.{idx}.mem.ptc")
        mem_pct_config = HardwareConfig(
            y_range=(0, 100), chart_name="GPU Memory Allocated (%)", chart_index=random_index()
        )
        # GPU 温度
        self.gpu_temp_key = generate_key("gpu.{idx}.temp")
        tem_config = HardwareConfig(chart_name="GPU Temperature (℃)", chart_index=random_index())
        # GPU 功耗
        self.gpu_power_key = generate_key("gpu.{idx}.power")
        power_config = HardwareConfig(chart_name="GPU Power Usage (W)", chart_index=random_index())
        # 每个GPU的配置信息
        self.per_gpu_configs = {
            self.gpu_mem_pct_key: [],
            self.gpu_temp_key: [],
            self.gpu_power_key: [],
            self.gpu_util_key: [],
        }
        self.handles = []
        for idx in range(count):
            metric_name = "GPU {idx}".format(idx=idx)
            self.per_gpu_configs[self.gpu_mem_pct_key].append(mem_pct_config.clone(metric_name=metric_name))
            self.per_gpu_configs[self.gpu_temp_key].append(tem_config.clone(metric_name=metric_name))
            self.per_gpu_configs[self.gpu_power_key].append(power_config.clone(metric_name=metric_name))
            self.per_gpu_configs[self.gpu_util_key].append(util_config.clone(metric_name=metric_name))

    @HardwareCollector.try_run()
    def get_gpu_config(self, key: str, idx: int) -> HardwareConfig:
        """
        获取 某个GPU的某个配置信息
        """
        return self.per_gpu_configs[key][idx]

    @HardwareCollector.try_run()
    def get_gpu_util(self, idx: int) -> HardwareInfo:
        """
        获取 GPU 利用率
        """
        handle = self.handles[idx]
        util_info = pynvml.nvmlDeviceGetUtilizationRates(handle)
        return {
            "key": self.gpu_util_key.format(idx=idx),
            "value": util_info.gpu,
            "name": "GPU {idx} Utilization (%)".format(idx=idx),
            "config": self.get_gpu_config(self.gpu_util_key, idx),
        }

    @HardwareCollector.try_run()
    def get_gpu_mem_pct(self, idx: int) -> HardwareInfo:
        """
        获取 GPU 内存使用率
        """
        handle = self.handles[idx]
        mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
        mem_pct = mem_info.used / mem_info.total * 100
        return {
            "key": self.gpu_mem_pct_key.format(idx=idx),
            "value": mem_pct,
            "name": "GPU {idx} Memory Allocated (%)".format(idx=idx),
            "config": self.get_gpu_config(self.gpu_mem_pct_key, idx),
        }

    @HardwareCollector.try_run()
    def get_gpu_temp(self, idx: int) -> HardwareInfo:
        """
        获取 GPU 温度
        """
        handle = self.handles[idx]
        temp_info = pynvml.nvmlDeviceGetTemperature(handle, pynvml.NVML_TEMPERATURE_GPU)
        return {
            "key": self.gpu_temp_key.format(idx=idx),
            "value": temp_info,
            "name": "GPU {idx} Temperature (℃)".format(idx=idx),
            "config": self.get_gpu_config(self.gpu_temp_key, idx),
        }

    @HardwareCollector.try_run()
    def get_gpu_power(self, idx: int) -> HardwareInfo:
        """
        获取 GPU 功耗
        """
        handle = self.handles[idx]
        # 功耗单位为mW，转换为W
        power_info = pynvml.nvmlDeviceGetPowerUsage(handle) / 1000
        return {
            "key": self.gpu_power_key.format(idx=idx),
            "value": power_info,
            "name": "GPU {idx} Power Usage (W)".format(idx=idx),
            "config": self.get_gpu_config(self.gpu_power_key, idx),
        }

    def collect(self) -> HardwareInfoList:
        """
        采集信息
        """
        result: HardwareInfoList = []
        for idx, handle in enumerate(self.handles):
            result.append(self.get_gpu_util(idx))
            result.append(self.get_gpu_mem_pct(idx))
            result.append(self.get_gpu_temp(idx))
            result.append(self.get_gpu_power(idx))
        return result

    def __del__(self):
        pynvml.nvmlShutdown()

    def before_collect_impl(self):
        # 低频采集下（30s以下），应该每次采集时都执行pynvml.nvmlInit()
        # 高频采集下（30s以上），应该在初始化时执行pynvml.nvmlInit()，在最后一次采集时执行pynvml.nvmlShutdown()
        # 在外部定时任务处，超过10次即变为低频采集，因此需要判断一下
        pynvml.nvmlInit()
        for i in range(pynvml.nvmlDeviceGetCount()):
            self.handles.append(pynvml.nvmlDeviceGetHandleByIndex(i))
        swanlog.debug("NVIDIA GPU nvml inited.")

    def after_collect_impl(self):
        pynvml.nvmlShutdown()
        self.handles.clear()
        swanlog.debug("NVIDIA GPU nvml shutdown.")
