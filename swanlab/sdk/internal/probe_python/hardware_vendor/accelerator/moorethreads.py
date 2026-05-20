"""
@author: cunyue
@file: moorethreads.py
@time: 2026/3/31 01:57
@description: Moore Threads GPU 信息采集模块

检测原理：
- 仅 Linux 支持，依赖摩尔线程 `mthreads-gmi` 的 JSON 输出。
- 静态信息使用 `mthreads-gmi -q --json`，顶层 `Driver Version` 为驱动版本，`GPU` 数组内每个对象包含
  `Index`、`Product Name`、`FB Memory Usage.Total`（`MiB` 或 `GiB`）。
- 利用率使用 `mthreads-gmi -q -d UTILIZATION --json`，解析 `GPU[].Utilization.Gpu`。
- 显存使用使用 `mthreads-gmi -q -d MEMORY --json`，解析 `GPU[].FB Memory Usage.Total/Used`（MiB）。
- 温度使用 `mthreads-gmi -q -d TEMPERATURE --json`，解析 `GPU[].Temperature.GPU Current Temp`。
- 功耗使用 `mthreads-gmi -q -d POWER --json`，解析 `GPU[].Power Readings["Power Draw "]`。
"""

import json
import math
import platform
import subprocess
from typing import Dict, List, Optional, Tuple

from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.typings.probe_python import (
    AcceleratorSnapshot,
    DeviceSnapshot,
    SystemScalar,
    SystemScalars,
    SystemShim,
)
from swanlab.sdk.typings.probe_python.hardware_vendor import AcceleratorProtocol, CollectResult
from swanlab.utils import generate_color


class MooreThreadsGPU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        mt_config = next(a for a in shim.accelerators if a.vendor == "moorethreads")
        self._indices = mt_config.device_indices
        driver, gpu_map = MooreThreadsGPU._map_mt_gpu()
        self._gpu_map = gpu_map
        max_mem_mb = 0
        for idx in self._indices:
            mem = int(gpu_map.get(str(idx), {}).get("memory", 0))
            if mem * 1024 > max_mem_mb:
                max_mem_mb = mem * 1024
        self._max_mem_mb = max_mem_mb

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize MooreThreads GPU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["MooreThreadsGPU", SystemScalars]]:
        mt_config = next((a for a in shim.accelerators if a.vendor == "moorethreads"), None)
        if mt_config is None:
            console.debug("No MooreThreads monitor config in shim")
            return None
        self = cls(shim)

        scalars: SystemScalars = []
        for color_idx, idx in enumerate(self._indices):
            color = generate_color(color_idx)

            scalars.append(
                SystemScalar(
                    key=f"gpu.{idx}.pct",
                    name=f"GPU {idx}",
                    chart_name="GPU Utilization (%)",
                    y_min=0,
                    y_max=100,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"gpu.{idx}.mem.pct",
                    name=f"GPU {idx}",
                    chart_name="GPU Memory Allocated (%)",
                    y_min=0,
                    y_max=100,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"gpu.{idx}.mem.value",
                    name=f"GPU {idx}",
                    chart_name="GPU Memory Allocated (MB)",
                    y_min=0,
                    y_max=self._max_mem_mb,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"gpu.{idx}.temp",
                    name=f"GPU {idx}",
                    chart_name="GPU Temperature (°C)",
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"gpu.{idx}.power",
                    name=f"GPU {idx}",
                    chart_name="GPU Power (W)",
                    color=color,
                )
            )

        console.debug(f"MooreThreads monitor initialized with {len(self._indices)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        results.extend(self._get_utilization())
        results.extend(self._get_memory())
        results.extend(self._get_temperature())
        results.extend(self._get_power())
        return results

    def _get_utilization(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.pct", math.nan))
        with safe.block(message="Failed to collect MooreThreads GPU utilization data", level="debug"):
            output = subprocess.run(
                ["mthreads-gmi", "-q", "-d", "UTILIZATION", "--json"],
                capture_output=True,
                text=True,
            ).stdout
            data = json.loads(output)
            for gpu_info in data.get("GPU", []):
                gpu_id = gpu_info["Index"]
                if gpu_id in self._indices:
                    gpu_util = gpu_info["Utilization"]["Gpu"]
                    if isinstance(gpu_util, str) and gpu_util.endswith("%"):
                        gpu_util = gpu_util[:-1]
                    idx_pos = self._indices.index(gpu_id)
                    results[idx_pos] = (f"gpu.{gpu_id}.pct", float(gpu_util))
        return results

    def _get_memory(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.mem.pct", math.nan))
            results.append((f"gpu.{idx}.mem.value", math.nan))
        with safe.block(message="Failed to collect MooreThreads GPU memory data", level="debug"):
            output = subprocess.run(
                ["mthreads-gmi", "-q", "-d", "MEMORY", "--json"],
                capture_output=True,
                text=True,
            ).stdout
            data = json.loads(output)
            for gpu_info in data.get("GPU", []):
                gpu_id = gpu_info["Index"]
                if gpu_id in self._indices:
                    mem_info = gpu_info["FB Memory Usage"]
                    mem_total = float(mem_info["Total"].replace("MiB", "").strip())
                    mem_used = float(mem_info["Used"].replace("MiB", "").strip())
                    pct = (mem_used / mem_total * 100) if mem_total > 0 else math.nan
                    idx_pos = self._indices.index(gpu_id)
                    results[idx_pos * 2] = (f"gpu.{gpu_id}.mem.pct", pct)
                    results[idx_pos * 2 + 1] = (f"gpu.{gpu_id}.mem.value", mem_used)
        return results

    def _get_temperature(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.temp", math.nan))
        with safe.block(message="Failed to collect MooreThreads GPU temperature data", level="debug"):
            output = subprocess.run(
                ["mthreads-gmi", "-q", "-d", "TEMPERATURE", "--json"],
                capture_output=True,
                text=True,
            ).stdout
            data = json.loads(output)
            for gpu_info in data.get("GPU", []):
                gpu_id = gpu_info["Index"]
                if gpu_id in self._indices:
                    temp_str = gpu_info["Temperature"]["GPU Current Temp"].replace("C", "").strip()
                    idx_pos = self._indices.index(gpu_id)
                    results[idx_pos] = (f"gpu.{gpu_id}.temp", float(temp_str))
        return results

    def _get_power(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.power", math.nan))
        with safe.block(message="Failed to collect MooreThreads GPU power data", level="debug"):
            output = subprocess.run(
                ["mthreads-gmi", "-q", "-d", "POWER", "--json"],
                capture_output=True,
                text=True,
            ).stdout
            data = json.loads(output)
            for gpu_info in data.get("GPU", []):
                gpu_id = gpu_info["Index"]
                if gpu_id in self._indices:
                    power_str = gpu_info["Power Readings"]["Power Draw "].strip().replace("W", "")
                    idx_pos = self._indices.index(gpu_id)
                    results[idx_pos] = (f"gpu.{gpu_id}.power", float(power_str))
        return results

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get MooreThreads GPU info")
    def get() -> Optional[AcceleratorSnapshot]:
        if platform.system() != "Linux":
            console.debug(
                f"MooreThreads GPU detection skipped: mthreads-gmi is only supported on Linux, current platform is {platform.system()}"
            )
            return None

        driver, gpu_map = MooreThreadsGPU._map_mt_gpu()
        if not gpu_map:
            console.debug("No MooreThreads GPU detected: mthreads-gmi -q --json returned no GPU devices")
            return None

        devices = []
        for gpu_id in sorted(gpu_map.keys(), key=int):
            info = gpu_map[gpu_id]
            devices.append(
                DeviceSnapshot(
                    index=int(gpu_id),
                    name=info.get("name", "MTT GPU"),
                    memory=int(info.get("memory", 0)),
                    memory_unit="GB",
                )
            )

        names = [d.name for d in devices]
        console.debug(f"Detected {len(devices)} MooreThreads GPU(s): {', '.join(names)}")
        return AcceleratorSnapshot(
            vendor="moorethreads",
            version=driver,
            devices=devices,
        )

    @staticmethod
    def _map_mt_gpu() -> Tuple[Optional[str], dict]:
        output_str = subprocess.run(["mthreads-gmi", "-q", "--json"], capture_output=True, check=True, text=True).stdout
        output_json = json.loads(output_str)
        driver = output_json.get("Driver Version")
        gpu_map: Dict[str, dict] = {}

        for gpu_info in output_json.get("GPU", []):
            gpu_id = str(gpu_info["Index"])
            gpu_name = gpu_info["Product Name"]
            gpu_memory = gpu_info["FB Memory Usage"]["Total"]
            if gpu_memory.endswith("MiB"):
                total = int(gpu_memory[:-3]) // 1024
            elif gpu_memory.endswith("GiB"):
                total = int(gpu_memory[:-3])
            else:
                total = 0
            gpu_map[gpu_id] = {"name": gpu_name, "memory": str(total)}

        return driver, gpu_map
