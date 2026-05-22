"""
@author: cunyue
@file: metax.py
@time: 2026/3/31 01:57
@description: MetaX GPU 信息采集模块

检测原理：
- 仅 Linux 支持，依赖沐曦 `mx-smi` 命令输出文本。
- 静态信息使用无参数 `mx-smi`：包含 `mx-smi` 的行末为驱动版本；包含 `MACA` 的行为 MACA 版本；
  包含 `MetaX` 且不含 `Management` 的行为设备型号；包含 `MiB` 的显存行形如 `used/total MiB`，
  解析 `/` 后的 total 并转换为 GB。
- 利用率使用 `mx-smi --show-usage`，解析同时包含 `GPU` 和 `%` 的行，倒数第二列为利用率。
- 显存使用使用 `mx-smi --show-memory`，解析 `vis_vram total` 与 `vis_vram used` 行。
- 温度使用 `mx-smi --show-temperature` 的 `hotspot` 行；功耗使用 `mx-smi --show-board-power`
  中同时包含 `Total` 和 `W` 的行，均取倒数第二列数值。
"""

import math
import platform
import subprocess
from typing import Dict, List, Optional, Tuple

from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.internal.probe_python.protocol import AcceleratorProtocol, CollectResult
from swanlab.sdk.internal.probe_python.typings import (
    AcceleratorSnapshot,
    DeviceSnapshot,
    SystemScalar,
    SystemScalars,
    SystemShim,
)
from swanlab.utils import generate_color


class MetaXGPU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        metax_config = next(a for a in shim.accelerators if a.vendor == "metax")
        self._indices = metax_config.device_indices
        driver, maca, gpu_map = MetaXGPU._map_metax_gpu()
        self._gpu_map = gpu_map
        max_mem_mb = 0
        for idx in self._indices:
            mem = gpu_map.get(idx, {}).get("memory", 0)
            if isinstance(mem, str):
                mem = int(mem)
            if mem * 1024 > max_mem_mb:
                max_mem_mb = mem * 1024
        self._max_mem_mb = max_mem_mb

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize MetaX GPU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["MetaXGPU", SystemScalars]]:
        metax_config = next((a for a in shim.accelerators if a.vendor == "metax"), None)
        if metax_config is None:
            console.debug("No MetaX monitor config in shim")
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
                    y_min=0,
                    color=color,
                )
            )

        console.debug(f"MetaX monitor initialized with {len(self._indices)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        with safe.block(message="Failed to collect MetaX GPU metrics", level="debug"):
            results.extend(self._get_utilization())
            results.extend(self._get_memory())
            results.extend(self._get_temperature())
            results.extend(self._get_power())
        return results

    def _get_utilization(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.pct", math.nan))
        with safe.block(message="Failed to collect MetaX GPU utilization data", level="debug"):
            output = subprocess.run(["mx-smi", "--show-usage"], capture_output=True, check=True, text=True).stdout
            index = 0
            for line in output.split("\n"):
                if "GPU" in line and "%" in line:
                    parts = line.split()
                    val = float(parts[-2]) if len(parts) >= 2 else math.nan
                    if index < len(results):
                        results[index] = (f"gpu.{self._indices[index]}.pct", val)
                    index += 1
        return results

    def _get_memory(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.mem.pct", math.nan))
            results.append((f"gpu.{idx}.mem.value", math.nan))
        with safe.block(message="Failed to collect MetaX GPU memory data", level="debug"):
            output = subprocess.run(["mx-smi", "--show-memory"], capture_output=True, text=True).stdout
            gpu_mem_total = None
            index = 0
            for line in output.split("\n"):
                if "vis_vram total" in line:
                    gpu_mem_total = float(line.split()[-2])
                if "vis_vram used" in line and gpu_mem_total is not None:
                    gpu_mem_used = float(line.split()[-2])
                    pct = gpu_mem_used / gpu_mem_total * 100 if gpu_mem_total > 0 else math.nan
                    mb = gpu_mem_used
                    if index < len(self._indices):
                        pos = index * 2
                        results[pos] = (f"gpu.{self._indices[index]}.mem.pct", pct)
                        results[pos + 1] = (f"gpu.{self._indices[index]}.mem.value", mb)
                    index += 1
                    gpu_mem_total = None
        return results

    def _get_temperature(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.temp", math.nan))
        with safe.block(message="Failed to collect MetaX GPU temperature data", level="debug"):
            output = subprocess.run(["mx-smi", "--show-temperature"], capture_output=True, check=True, text=True).stdout
            index = 0
            for line in output.split("\n"):
                if "hotspot" in line:
                    parts = line.split()
                    val = float(parts[-2]) if len(parts) >= 2 else math.nan
                    if index < len(results):
                        results[index] = (f"gpu.{self._indices[index]}.temp", val)
                    index += 1
        return results

    def _get_power(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.power", math.nan))
        with safe.block(message="Failed to collect MetaX GPU power data", level="debug"):
            output = subprocess.run(["mx-smi", "--show-board-power"], capture_output=True, check=True, text=True).stdout
            index = 0
            for line in output.split("\n"):
                if "Total" in line and "W" in line:
                    parts = line.split()
                    val = float(parts[-2]) if len(parts) >= 2 else math.nan
                    if index < len(results):
                        results[index] = (f"gpu.{self._indices[index]}.power", val)
                    index += 1
        return results

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get MetaX GPU info")
    def get() -> Optional[AcceleratorSnapshot]:
        if platform.system() != "Linux":
            console.debug(
                f"MetaX GPU detection skipped: mx-smi is only supported on Linux, current platform is {platform.system()}"
            )
            return None

        driver, maca_version, gpu_map = MetaXGPU._map_metax_gpu()
        if not gpu_map:
            console.debug("No MetaX GPU detected: mx-smi output did not include MetaX device rows")
            return None

        devices = []
        for gpu_id in sorted(gpu_map.keys()):
            info = gpu_map[gpu_id]
            devices.append(
                DeviceSnapshot(
                    index=int(gpu_id),
                    name=info.get("name", "MetaX GPU"),
                    memory=int(info.get("memory", 0)),
                    memory_unit="GB",
                )
            )

        names = [d.name for d in devices]
        console.debug(f"Detected {len(devices)} MetaX GPU(s): {', '.join(names)}")
        return AcceleratorSnapshot(
            vendor="metax",
            version=driver,
            devices=devices,
        )

    @staticmethod
    def _map_metax_gpu() -> Tuple[Optional[str], Optional[str], Dict[int, dict]]:
        output_str = subprocess.run(["mx-smi"], capture_output=True, check=True, text=True).stdout

        driver = None
        maca_version = None
        gpu_map: Dict[int, dict] = {}
        index = 0

        for line in output_str.split("\n"):
            if "mx-smi" in line:
                driver = line.split()[-1]
            if "MACA" in line:
                maca_version = line.split()[3]
            elif "MetaX" in line and "Management" not in line:
                parts = line.split()
                if len(parts) >= 4:
                    gpu_map[index] = {"name": parts[3]}
            elif "MiB" in line and index in gpu_map:
                with safe.block(message="Failed to parse MetaX GPU memory info", level="debug"):
                    parts = line.split()
                    if len(parts) >= 4:
                        gpu_map[index]["memory"] = int(parts[-4].split("/")[-1]) // 1024
                        index += 1

        return driver, maca_version, gpu_map
