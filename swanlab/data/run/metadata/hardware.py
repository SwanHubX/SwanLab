"""
@author: cunyue
@file: hardware.py
@time: 2024/11/18 15:12
@description: 硬件信息采集
"""

import json
import multiprocessing
import platform
import subprocess

import psutil
import pynvml


def get_hardware_info():
    """
    采集硬件信息，包括CPU、GPU、内存、硬盘等
    """
    info = {
        "memory": get_memory_size(),
        "cpu": get_cpu_info(),
        "gpu": {
            "nvidia": get_nvidia_gpu_info(),
        },
        "soc": {
            "apple": get_apple_chip_info(),
        },
    }
    return info


# ---------------------------------- cpu信息 ----------------------------------


def get_cpu_info():
    """获取 CPU 信息"""
    info = {"brand": None, "cores": None}

    # 获取 CPU 品牌, 根据不同操作系统调用不同的函数
    if platform.system() == "Windows":
        info["brand"] = get_cpu_brand_windows()
    elif platform.system() == "Linux":
        info["brand"] = get_cpu_brand_linux()
    else:
        # 其他情况，暂时不支持
        # 苹果芯片单独处理
        return None
    try:
        # 获取 CPU 核心数
        info["cores"] = multiprocessing.cpu_count()
    except Exception:  # noqa
        pass

    return info


def get_cpu_brand_windows():
    try:
        # 使用 WMIC 命令获取 CPU 品牌
        result = subprocess.run(["wmic", "cpu", "get", "name"], capture_output=True, text=True)
        cpu_brand = result.stdout.strip().split("\n")[-1].strip()
        return cpu_brand
    except Exception:  # noqa
        return None


def get_cpu_brand_linux():
    try:
        # 使用 lscpu 命令获取 CPU 品牌
        result = subprocess.run(["lscpu"], capture_output=True, text=True)
        for line in result.stdout.split("\n"):
            if "Model name:" in line:
                cpu_brand = line.split(":")[1].strip()
                return cpu_brand
        return None
    except Exception:  # noqa
        return None


# ---------------------------------- 内存信息 ----------------------------------


def get_memory_size():
    """获取内存大小"""
    try:
        # 获取系统总内存大小
        mem = psutil.virtual_memory()
        total_memory = round(mem.total / (1024**3))  # 单位为GB
        return total_memory
    except Exception:  # noqa
        return


# ---------------------------------- gpu信息 ----------------------------------


def get_nvidia_gpu_info():
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
        return None

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
        return info


# ---------------------------------- apple信息 ----------------------------------


def get_apple_chip_info():
    if "mac" not in platform.platform().lower():
        return None
    info = {"cpu": None, "gpu": None, "memory": None, "type": None}

    # 使用system_profiler命令以JSON格式获取GPU信息
    try:
        result = subprocess.run(["system_profiler", "SPHardwareDataType", "-json"], capture_output=True, text=True)
        gpu_name = json.loads(result.stdout)["SPHardwareDataType"][0]["chip_type"]
        memory = json.loads(result.stdout)["SPHardwareDataType"][0]["physical_memory"]
        memory = str(memory).lower().replace("gb", "")
        # TODO: 获取GPU信息
        info["type"] = gpu_name
        info["memory"] = memory
    except Exception:  # noqa
        return None
    try:
        info["cpu"] = multiprocessing.cpu_count()
    except Exception:  # noqa
        pass
    return info
