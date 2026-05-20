"""
@author: cunyue
@file: hygon.py
@time: 2026/3/31 01:58
@description: 海光 Hygon DCU 信息采集模块

检测原理：
- 仅 Linux 支持，依赖海光 `hy-smi` 命令及其 JSON 输出。
- 驱动版本使用 `hy-smi --showdriverversion`，解析包含 `Driver Version` 的行。
- 静态设备信息使用 `hy-smi --showproductname --json` 的 `cardX -> Card Series`，以及
  `hy-smi --showmemavailable --json` 的 `cardX -> Available memory size (MiB)`；设备索引按 JSON 遍历顺序映射为 0,1,2...
- 利用率使用 `hy-smi --showuse --json`，解析 `DCU use (%)`。
- 显存使用率使用 `hy-smi --showmemuse --json`，解析 `DCU memory use (%)`；显存使用量按比例乘以最大显存 MB。
- 温度使用 `hy-smi --showtemp --json`，解析 `Temperature (Sensor junction) (C)`；功耗使用
  `hy-smi --showpower --json`，解析 `Average Graphics Package Power (W)`。
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


class HygonDCU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        hygon_config = next(a for a in shim.accelerators if a.vendor == "hygon")
        self._indices = hygon_config.device_indices
        driver, dcu_map = HygonDCU._map_hygon_dcu()
        self._dcu_map = dcu_map
        max_mem_mb = 0
        for idx in self._indices:
            dcu_id = str(idx)
            if dcu_id in dcu_map:
                mem_str = dcu_map[dcu_id].get("memory", "0GB")
                digits = [c for c in mem_str if c.isdigit()]
                mem_val = int("".join(digits)) if digits else 0
                if mem_val * 1024 > max_mem_mb:
                    max_mem_mb = mem_val * 1024
        self._max_mem_mb = max_mem_mb

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize Hygon DCU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["HygonDCU", SystemScalars]]:
        hygon_config = next((a for a in shim.accelerators if a.vendor == "hygon"), None)
        if hygon_config is None:
            console.debug("No Hygon DCU monitor config in shim")
            return None
        self = cls(shim)

        scalars: SystemScalars = []
        for color_idx, idx in enumerate(self._indices):
            color = generate_color(color_idx)

            scalars.append(
                SystemScalar(
                    key=f"dcu.{idx}.pct",
                    name=f"DCU {idx}",
                    chart_name="DCU Utilization (%)",
                    y_min=0,
                    y_max=100,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"dcu.{idx}.mem.pct",
                    name=f"DCU {idx}",
                    chart_name="DCU Memory Allocated (%)",
                    y_min=0,
                    y_max=100,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"dcu.{idx}.mem.value",
                    name=f"DCU {idx}",
                    chart_name="DCU Memory Allocated (MB)",
                    y_min=0,
                    y_max=self._max_mem_mb,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"dcu.{idx}.temp",
                    name=f"DCU {idx}",
                    chart_name="DCU Temperature (°C)",
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"dcu.{idx}.power",
                    name=f"DCU {idx}",
                    chart_name="DCU Power (W)",
                    color=color,
                )
            )

        console.debug(f"Hygon DCU monitor initialized with {len(self._indices)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        with safe.block(message="Failed to collect Hygon DCU metrics", level="debug"):
            results.extend(self._get_utilization())
            results.extend(self._get_memory())
            results.extend(self._get_temperature())
            results.extend(self._get_power())
        return results

    def _get_utilization(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"dcu.{idx}.pct", math.nan))
        with safe.block(message="Failed to collect Hygon DCU utilization data", level="debug"):
            output = subprocess.run(["hy-smi", "--showuse", "--json"], capture_output=True, text=True).stdout
            data = json.loads(output)
            for idx, (dcu_key, dcu_info) in enumerate(data.items()):
                if idx in self._indices:
                    results[self._indices.index(idx)] = (
                        f"dcu.{idx}.pct",
                        float(dcu_info["DCU use (%)"]),
                    )
        return results

    def _get_memory(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"dcu.{idx}.mem.pct", math.nan))
            results.append((f"dcu.{idx}.mem.value", math.nan))
        with safe.block(message="Failed to collect Hygon DCU memory data", level="debug"):
            output = subprocess.run(["hy-smi", "--showmemuse", "--json"], capture_output=True, text=True).stdout
            data = json.loads(output)
            for idx, (dcu_key, dcu_info) in enumerate(data.items()):
                if idx in self._indices:
                    rate = float(dcu_info["DCU memory use (%)"])
                    mb = rate * 0.01 * self._max_mem_mb
                    idx_pos = self._indices.index(idx)
                    results[idx_pos * 2] = (f"dcu.{idx}.mem.pct", rate)
                    results[idx_pos * 2 + 1] = (f"dcu.{idx}.mem.value", mb)
        return results

    def _get_temperature(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"dcu.{idx}.temp", math.nan))
        with safe.block(message="Failed to collect Hygon DCU temperature data", level="debug"):
            output = subprocess.run(["hy-smi", "--showtemp", "--json"], capture_output=True, text=True).stdout
            data = json.loads(output)
            for idx, (dcu_key, dcu_info) in enumerate(data.items()):
                if idx in self._indices:
                    results[self._indices.index(idx)] = (
                        f"dcu.{idx}.temp",
                        float(dcu_info["Temperature (Sensor junction) (C)"]),
                    )
        return results

    def _get_power(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            results.append((f"dcu.{idx}.power", math.nan))
        with safe.block(message="Failed to collect Hygon DCU power data", level="debug"):
            output = subprocess.run(["hy-smi", "--showpower", "--json"], capture_output=True, text=True).stdout
            data = json.loads(output)
            for idx, (dcu_key, dcu_info) in enumerate(data.items()):
                if idx in self._indices:
                    results[self._indices.index(idx)] = (
                        f"dcu.{idx}.power",
                        float(dcu_info["Average Graphics Package Power (W)"]),
                    )
        return results

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get Hygon DCU info")
    def get() -> Optional[AcceleratorSnapshot]:
        if platform.system() != "Linux":
            console.debug(
                f"Hygon DCU detection skipped: hy-smi is only supported on Linux, current platform is {platform.system()}"
            )
            return None

        driver, dcu_map = HygonDCU._map_hygon_dcu()
        if not dcu_map:
            console.debug("No Hygon DCU detected: hy-smi product/memory queries returned no DCU devices")
            return None

        devices = []
        for idx, dcu_id in enumerate(sorted(dcu_map.keys(), key=int)):
            info = dcu_map[dcu_id]
            mem_str = info.get("memory", "0GB")
            digits = [c for c in mem_str if c.isdigit()]
            mem_val = int("".join(digits)) if digits else 0
            devices.append(
                DeviceSnapshot(
                    index=int(dcu_id),
                    name=info.get("name", "Hygon DCU"),
                    memory=mem_val,
                    memory_unit="GB",
                )
            )

        names = [d.name for d in devices]
        console.debug(f"Detected {len(devices)} Hygon DCU device(s): {', '.join(names)}")
        return AcceleratorSnapshot(
            vendor="hygon",
            version=driver,
            devices=devices,
        )

    @staticmethod
    def _map_hygon_dcu() -> Tuple[Optional[str], dict]:
        driver_version = None
        with safe.block(message="Failed to get Hygon DCU driver version", level="debug"):
            driver_info = subprocess.run(
                ["hy-smi", "--showdriverversion"], capture_output=True, check=True, text=True
            ).stdout
            for line in driver_info.split("\n"):
                if "Driver Version" in line:
                    driver_version = line.split(":")[-1].strip()
                    break

        dcu_map: Dict[str, dict] = {}
        with safe.block(message="Failed to map Hygon DCU devices", level="debug"):
            product_map = json.loads(
                subprocess.run(
                    ["hy-smi", "--showproductname", "--json"], capture_output=True, check=True, text=True
                ).stdout
            )
            available_mem_map = json.loads(
                subprocess.run(
                    ["hy-smi", "--showmemavailable", "--json"], capture_output=True, check=True, text=True
                ).stdout
            )

            for idx, device_id in enumerate(product_map.keys()):
                product_info = product_map.get(device_id, {})
                device_name = product_info.get("Card Series", "Unknown")
                available_info = available_mem_map.get(device_id, {})
                available_mem = available_info.get("Available memory size (MiB)", "0")
                total_mem_gb = round(int(available_mem) / 1024)
                dcu_map[str(idx)] = {"name": device_name, "memory": f"{total_mem_gb}GB"}

        return driver_version, dcu_map
