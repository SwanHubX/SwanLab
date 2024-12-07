"""
@author: cunyue
@file: nvidia.py
@time: 2024/12/3 20:12
@description: NVIDIA GPU信息采集
"""

import subprocess

import pynvml

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult


def get_nvidia_gpu_info() -> HardwareFuncResult:
    """获取 GPU 信息"""

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

    except pynvml.NVMLError:
        pass
    finally:
        # 结束 NVML
        pynvml.nvmlShutdown()
        return info, None
