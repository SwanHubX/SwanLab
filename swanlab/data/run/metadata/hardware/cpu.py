"""
@author: cunyue
@file: cpu.py
@time: 2024/12/3 20:12
@description: CPU信息采集
"""

import multiprocessing
import platform
import subprocess

from swanlab.data.run.metadata.hardware.type import HardwareFuncResult


def get_cpu_info() -> HardwareFuncResult:
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
        return None, []
    try:
        # 获取 CPU 核心数
        info["cores"] = multiprocessing.cpu_count()
    except Exception:  # noqa
        pass

    return info, []


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
