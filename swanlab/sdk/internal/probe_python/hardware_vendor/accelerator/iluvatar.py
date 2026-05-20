"""
@author: cunyue
@file: iluvatar.py
@time: 2026/3/31 01:58
@description: Iluvatar GPU 信息采集模块

检测原理：
- 仅 Linux 支持，依赖天数智芯 `ixsmi` 命令的 CSV 查询输出。
- 静态信息使用 `ixsmi --query-gpu=index,driver_version,name,memory.total --format=csv`，
  第一行为表头，后续每行格式为 `index, driver_version, name, memory.total`。
- `memory.total` 可能以 `MiB` 或 `GiB` 结尾，统一转换为 GB。
- 动态指标分别使用 `--query-gpu=utilization.gpu|utilization.memory|memory.used|temperature.gpu|board.power.draw`，
  每个 CSV 输出跳过表头，解析第一列并去掉 `%`、`MiB`、`C`、`W` 单位。
"""

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


class IluvatarGPU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        iluvatar_config = next(a for a in shim.accelerators if a.vendor == "iluvatar")
        self._indices = iluvatar_config.device_indices
        _, gpu_map = IluvatarGPU._map_iluvatar_gpu()
        self._gpu_map = gpu_map
        max_mem_mb = 0
        for idx in self._indices:
            gpu_id = str(idx)
            if gpu_id in gpu_map:
                mem = int(gpu_map[gpu_id].get("memory", 0))
                if mem * 1024 > max_mem_mb:
                    max_mem_mb = mem * 1024
        self._max_mem_mb = max_mem_mb

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize Iluvatar GPU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["IluvatarGPU", SystemScalars]]:
        iluvatar_config = next((a for a in shim.accelerators if a.vendor == "iluvatar"), None)
        if iluvatar_config is None:
            console.debug("No Iluvatar monitor config in shim")
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

        console.debug(f"Iluvatar monitor initialized with {len(self._indices)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        with safe.block(message="Failed to collect Iluvatar GPU metrics", level="debug"):
            results.extend(self._collect_utilization())
            results.extend(self._collect_memory())
            results.extend(self._collect_temperature())
            results.extend(self._collect_power())
        return results

    def _collect_utilization(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.pct", math.nan))
        with safe.block(message="Failed to collect Iluvatar GPU utilization", level="debug"):
            gpu_utils = self._run_ixsmi_query("utilization.gpu")
            for i, val in enumerate(gpu_utils):
                if i < len(self._indices):
                    results[i] = (f"gpu.{self._indices[i]}.pct", val)
        return results

    def _collect_memory(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.mem.pct", math.nan))
            results.append((f"gpu.{idx}.mem.value", math.nan))
        with safe.block(message="Failed to collect Iluvatar GPU memory", level="debug"):
            mem_rates = self._run_ixsmi_query("utilization.memory")
            mem_used = self._run_ixsmi_query("memory.used")
            for i in range(min(len(self._indices), len(mem_rates), len(mem_used))):
                results[i * 2] = (f"gpu.{self._indices[i]}.mem.pct", mem_rates[i])
                results[i * 2 + 1] = (f"gpu.{self._indices[i]}.mem.value", mem_used[i])
        return results

    def _collect_temperature(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.temp", math.nan))
        with safe.block(message="Failed to collect Iluvatar GPU temperature", level="debug"):
            temps = self._run_ixsmi_query("temperature.gpu")
            for i, val in enumerate(temps):
                if i < len(self._indices):
                    results[i] = (f"gpu.{self._indices[i]}.temp", val)
        return results

    def _collect_power(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"gpu.{idx}.power", math.nan))
        with safe.block(message="Failed to collect Iluvatar GPU power", level="debug"):
            powers = self._run_ixsmi_query("board.power.draw")
            for i, val in enumerate(powers):
                if i < len(self._indices):
                    results[i] = (f"gpu.{self._indices[i]}.power", val)
        return results

    @staticmethod
    def _run_ixsmi_query(query_type: str) -> List[float]:
        output = subprocess.run(
            ["ixsmi", f"--query-gpu={query_type}", "--format=csv"],
            capture_output=True,
            check=True,
            text=True,
            encoding="utf-8",
        ).stdout
        lines = output.strip().split("\n")
        results: List[float] = []
        for line in lines[1:]:
            value_str = line.strip().split(",")[0].strip()
            if query_type == "utilization.gpu" or query_type == "utilization.memory":
                value = float(value_str.replace("%", "").strip())
            elif query_type == "temperature.gpu":
                value = float(value_str.replace("C", "").strip())
            elif query_type == "memory.used":
                value = float(value_str.replace("MiB", "").strip())
            elif query_type == "board.power.draw":
                value = float(value_str.replace("W", "").strip())
            else:
                value = 0.0
            results.append(value)
        return results

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get Iluvatar GPU info")
    def get() -> Optional[AcceleratorSnapshot]:
        if platform.system() != "Linux":
            console.debug(
                f"Iluvatar GPU detection skipped: ixsmi is only supported on Linux, current platform is {platform.system()}"
            )
            return None

        driver, gpu_map = IluvatarGPU._map_iluvatar_gpu()
        if not gpu_map:
            console.debug("No Iluvatar GPU detected: ixsmi query returned no GPU rows")
            return None

        devices = []
        for gpu_id in sorted(gpu_map.keys(), key=int):
            info = gpu_map[gpu_id]
            devices.append(
                DeviceSnapshot(
                    index=int(gpu_id),
                    name=info.get("name", "Iluvatar GPU"),
                    memory=int(info.get("memory", 0)),
                    memory_unit="GB",
                )
            )

        names = [d.name for d in devices]
        console.debug(f"Detected {len(devices)} Iluvatar GPU(s): {', '.join(names)}")
        return AcceleratorSnapshot(
            vendor="iluvatar",
            version=driver,
            devices=devices,
        )

    @staticmethod
    def _map_iluvatar_gpu() -> Tuple[Optional[str], dict]:
        output_str = None
        with safe.block(message="Failed to query Iluvatar GPU info via ixsmi", level="debug"):
            output_str = subprocess.run(
                ["ixsmi", "--query-gpu=index,driver_version,name,memory.total", "--format=csv"],
                capture_output=True,
                check=True,
                text=True,
                encoding="utf-8",
            ).stdout

        if output_str is None:
            return None, {}

        lines = output_str.strip().split("\n")
        if len(lines) < 2:
            return None, {}

        driver = None
        gpu_map: Dict[str, dict] = {}

        for line in lines[1:]:
            parts = [part.strip() for part in line.split(",")]
            if len(parts) >= 4:
                gpu_id = parts[0]
                if driver is None:
                    driver = parts[1]
                gpu_name = parts[2]
                gpu_memory_str = parts[3]
                if "MiB" in gpu_memory_str:
                    memory_gb = str(int(gpu_memory_str.replace(" MiB", "").strip()) // 1024)
                elif "GiB" in gpu_memory_str:
                    memory_gb = str(int(gpu_memory_str.replace(" GiB", "").strip()))
                else:
                    memory_gb = gpu_memory_str
                gpu_map[gpu_id] = {"name": gpu_name, "memory": memory_gb}

        return driver, gpu_map
