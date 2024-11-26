"""
@author: cunyue
@file: hardware.py
@time: 2024/11/18 15:12
@description: 硬件信息采集
"""

import json
import multiprocessing
import os
import platform
import subprocess

import psutil
import pynvml

from .utils import ascend


def get_hardware_info():
    """
    采集硬件信息，包括CPU、GPU、内存、硬盘等
    """
    gpu = filter_none({"nvidia": get_nvidia_gpu_info()})
    npu = filter_none({"ascend": get_ascend_npu_info()})
    soc = filter_none({"apple": get_apple_chip_info()})

    info = {
        "memory": get_memory_size(),
        "cpu": get_cpu_info(),
        "gpu": gpu,
        "npu": npu,
        "soc": soc,
    }
    return filter_none(info, fallback={})


def filter_none(data, fallback=None):
    """
    过滤掉字典中值为None的键值对，只对字典有效
    """
    if isinstance(data, dict):
        data = {k: v for k, v in data.items() if v is not None}
        if all(v is None for v in data.values()):
            return fallback
    return data


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


# ---------------------------------- npu信息 ----------------------------------


def get_ascend_npu_info():
    """
    获取华为昇腾NPU信息，包括驱动版本、设备信息等
    目前的信息统计粒度只到NPU ID级别，没有到Chip ID级别
    """
    # ascend芯片只支持Linux系统
    if platform.system() != "Linux":
        return None
    # /dev目录下没有davinci*设备文件，跳过
    # 其实理论上davinci后接数字，代表此设备id，但是官方文档也没明确写，以防万一还是不这么干了
    if not list(filter(lambda x: x.startswith("davinci"), os.listdir("/dev"))):
        return None
    info = {"driver": None, "npu": None}
    try:
        # 获取NPU驱动版本
        info["driver"] = ascend.get_version()
        # 获取所有NPU设备ID
        npu_map = ascend.map_npu()
        for npu_id in npu_map:
            for chip_id in npu_map[npu_id]:
                chip_info = npu_map[npu_id][chip_id]
                usage = ascend.get_chip_usage(npu_id, chip_id)
                if info["npu"] is None:
                    info["npu"] = {}
                if npu_id not in info["npu"]:
                    info["npu"][npu_id] = {}
                info["npu"][npu_id][chip_id] = {**chip_info, "usage": usage}
    except Exception:  # noqa
        if all(v is None for v in info.values()):
            return None
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
