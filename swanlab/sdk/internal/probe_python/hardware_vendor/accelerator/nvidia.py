"""
@author: cunyue
@file: nvidia.py
@time: 2026/3/31 01:58
@description: NVIDIA GPU 信息采集模块

检测原理：
- 首选 NVIDIA 官方 NVML（Python 包 pynvml）查询设备，NVML 是 nvidia-smi 的底层管理接口之一。
- 静态信息通过 nvmlSystemGetDriverVersion、nvmlDeviceGetCount、nvmlDeviceGetName、
  nvmlDeviceGetMemoryInfo 采集。
- CUDA 版本补充查询先使用 `nvcc --version`，解析形如 `Cuda compilation tools, release 12.4, V...` 的
  release 字段；失败后使用 `nvidia-smi` 表头中的 `CUDA Version: 12.4` 字段兜底。
- 动态指标来自 NVML：nvmlDeviceGetUtilizationRates(handle).gpu/memory、
  nvmlDeviceGetMemoryInfo(handle).used/total、nvmlDeviceGetTemperature、nvmlDeviceGetPowerUsage。
"""

import subprocess
from typing import Optional, Tuple

import pynvml

from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.typings.probe_python import (
    AcceleratorSnapshot,
    DeviceSnapshot,
    SystemScalar,
    SystemScalars,
    SystemShim,
)
from swanlab.sdk.typings.probe_python.hardware_vendor import AcceleratorProtocol
from swanlab.utils import generate_color


class NvidiaGPU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        pynvml.nvmlInit()
        nvidia_config = next(a for a in shim.accelerators if a.vendor == "nvidia")
        self._indices = nvidia_config.device_indices
        self._handles = [pynvml.nvmlDeviceGetHandleByIndex(i) for i in self._indices]
        max_mem_mb = 0
        for h in self._handles:
            total = int(pynvml.nvmlDeviceGetMemoryInfo(h).total) >> 20
            if total > max_mem_mb:
                max_mem_mb = total
        self._max_mem_mb = max_mem_mb

    def __del__(self):
        try:
            pynvml.nvmlShutdown()
        except Exception:
            pass

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize NVIDIA GPU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["NvidiaGPU", SystemScalars]]:
        nvidia_config = next((a for a in shim.accelerators if a.vendor == "nvidia"), None)
        if nvidia_config is None:
            console.debug("No NVIDIA monitor config in shim")
            return None
        try:
            self = cls(shim)
        except pynvml.NVMLError:
            console.debug("NVIDIA monitor initialization skipped: pynvml.nvmlInit() failed")
            return None

        scalars: SystemScalars = []
        for color_idx, idx in enumerate(self._indices):
            handle = self._handles[color_idx]
            color = generate_color(color_idx)

            util = SystemScalar(
                key=f"gpu.{idx}.pct",
                name=f"GPU {idx}",
                chart_name="GPU Utilization (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(util)
            self._handlers.append(
                (f"gpu.{idx}.pct", lambda h=handle: float(pynvml.nvmlDeviceGetUtilizationRates(h).gpu))
            )

            mem_pct = SystemScalar(
                key=f"gpu.{idx}.mem.pct",
                name=f"GPU {idx}",
                chart_name="GPU Memory Allocated (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(mem_pct)
            self._handlers.append(
                (
                    f"gpu.{idx}.mem.pct",
                    lambda h=handle: (
                        float(pynvml.nvmlDeviceGetMemoryInfo(h).used)
                        / float(pynvml.nvmlDeviceGetMemoryInfo(h).total)
                        * 100
                    ),
                )
            )

            mem_value = SystemScalar(
                key=f"gpu.{idx}.mem.value",
                name=f"GPU {idx}",
                chart_name="GPU Memory Allocated (MB)",
                y_min=0,
                y_max=self._max_mem_mb,
                color=color,
            )
            scalars.append(mem_value)
            self._handlers.append(
                (f"gpu.{idx}.mem.value", lambda h=handle: int(pynvml.nvmlDeviceGetMemoryInfo(h).used) >> 20)
            )

            temp = SystemScalar(
                key=f"gpu.{idx}.temp",
                name=f"GPU {idx}",
                chart_name="GPU Temperature (°C)",
                color=color,
            )
            scalars.append(temp)
            self._handlers.append(
                (
                    f"gpu.{idx}.temp",
                    lambda h=handle: int(pynvml.nvmlDeviceGetTemperature(h, pynvml.NVML_TEMPERATURE_GPU)),
                )
            )

            power = SystemScalar(
                key=f"gpu.{idx}.power",
                name=f"GPU {idx}",
                chart_name="GPU Power Usage (W)",
                color=color,
            )
            scalars.append(power)
            self._handlers.append(
                (f"gpu.{idx}.power", lambda h=handle: float(pynvml.nvmlDeviceGetPowerUsage(h)) / 1000)
            )

            mem_time = SystemScalar(
                key=f"gpu.{idx}.mem.time",
                name=f"GPU {idx}",
                chart_name="GPU Time Spent Accessing Memory (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(mem_time)
            self._handlers.append(
                (f"gpu.{idx}.mem.time", lambda h=handle: float(pynvml.nvmlDeviceGetUtilizationRates(h).memory))
            )

        console.debug(f"NVIDIA monitor initialized with {len(self._indices)} device(s)")
        return self, scalars

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get NVIDIA GPU info")
    def get() -> Optional[AcceleratorSnapshot]:
        try:
            pynvml.nvmlInit()
        except pynvml.NVMLError:
            console.debug(
                "NVIDIA GPU detection skipped: pynvml.nvmlInit() failed; NVIDIA driver/NVML may be unavailable"
            )
            return None
        try:
            driver = pynvml.nvmlSystemGetDriverVersion()
            if isinstance(driver, bytes):
                driver = driver.decode("utf-8")
            count = pynvml.nvmlDeviceGetCount()
            if count == 0:
                console.debug("No NVIDIA GPU detected: NVML initialized but reported 0 devices")
            cuda_version = NvidiaGPU._get_cuda_version()
            devices = []
            for i in range(count):
                handle = pynvml.nvmlDeviceGetHandleByIndex(i)
                name = pynvml.nvmlDeviceGetName(handle)
                if isinstance(name, bytes):
                    name = name.decode("utf-8")
                mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                mem_gb = round(int(mem_info.total) / 1024**3)
                devices.append(DeviceSnapshot(index=i, name=name, memory=mem_gb, memory_unit="GB"))
            names = [d.name for d in devices]
            console.debug(f"Detected {count} NVIDIA GPU(s): {', '.join(names)}")
            return AcceleratorSnapshot(
                vendor="nvidia",
                version=driver,
                cuda_version=cuda_version,
                devices=devices,
            )
        finally:
            pynvml.nvmlShutdown()

    @staticmethod
    def _get_cuda_version() -> Optional[str]:
        return NvidiaGPU._get_cuda_version_from_nvcc() or NvidiaGPU._get_cuda_version_from_smi()

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get CUDA version via nvcc")
    def _get_cuda_version_from_nvcc() -> Optional[str]:
        output = subprocess.check_output(["nvcc", "--version"]).decode("utf-8")
        for line in output.split("\n"):
            if "release" in line:
                return line.split("release")[-1].strip().split(" ")[0][:-1]
        return None

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get CUDA version via nvidia-smi")
    def _get_cuda_version_from_smi() -> Optional[str]:
        cmd = "nvidia-smi | grep 'CUDA Version' | awk -F'CUDA Version: ' '{print $2}' | awk '{print $1}'"
        output = subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
        return output or None
