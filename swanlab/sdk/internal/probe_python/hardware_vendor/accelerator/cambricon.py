"""
@author: cunyue
@file: cambricon.py
@time: 2026/3/31 01:58
@description: 寒武纪 Cambricon MLU 信息采集模块

检测原理：
- 仅 Linux 支持，依赖寒武纪 `cnmon` 命令行工具输出文本。
- 静态信息使用 `cnmon info -m`，按段落解析：`Card ... <id>` 开启一个设备段，
  `Product ... <name>` 为型号，`Driver ... <version>` 为驱动版本，`Total ... <MiB>` 为总显存。
- 利用率使用 `cnmon info -u`，解析包含 `mlu average` 的行，冒号后的百分比为 MLU 使用率。
- 显存使用使用 `cnmon info -m` 的 `physical memory usage` 段，其后 `used: <MiB>` 行为已用显存。
- 温度使用 `cnmon info -e` 的 `temperature` 段，解析 `chip: <C>`；功耗使用 `cnmon info -p`
  的 `power` 段，解析 `usage: <W>`。
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


class CambriconMLU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        cambricon_config = next(a for a in shim.accelerators if a.vendor == "cambricon")
        self._indices = cambricon_config.device_indices
        self._mlu_info = CambriconMLU._map_mlu_raw()[1]

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize Cambricon MLU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["CambriconMLU", SystemScalars]]:
        cambricon_config = next((a for a in shim.accelerators if a.vendor == "cambricon"), None)
        if cambricon_config is None:
            console.debug("No Cambricon MLU monitor config in shim")
            return None
        self = cls(shim)

        max_mem_mb = 0
        for idx in cambricon_config.device_indices:
            mlu_id = str(idx)
            if mlu_id in self._mlu_info:
                mem = int(self._mlu_info[mlu_id].get("memory", 0))
                if mem * 1024 > max_mem_mb:
                    max_mem_mb = mem * 1024

        scalars: SystemScalars = []
        for color_idx, idx in enumerate(cambricon_config.device_indices):
            color = generate_color(color_idx)

            util = SystemScalar(
                key=f"mlu.{idx}.pct",
                name=f"MLU {idx}",
                chart_name="MLU Utilization (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(util)

            mem_pct = SystemScalar(
                key=f"mlu.{idx}.mem.pct",
                name=f"MLU {idx}",
                chart_name="MLU Memory Allocated (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(mem_pct)

            mem_value = SystemScalar(
                key=f"mlu.{idx}.mem.value",
                name=f"MLU {idx}",
                chart_name="MLU Memory Allocated (MB)",
                y_min=0,
                y_max=max_mem_mb,
                color=color,
            )
            scalars.append(mem_value)

            temp = SystemScalar(
                key=f"mlu.{idx}.temp",
                name=f"MLU {idx}",
                chart_name="MLU Temperature (°C)",
                color=color,
            )
            scalars.append(temp)

            power = SystemScalar(
                key=f"mlu.{idx}.power",
                name=f"MLU {idx}",
                chart_name="MLU Power (W)",
                y_min=0,
                color=color,
            )
            scalars.append(power)

        console.debug(f"Cambricon MLU monitor initialized with {len(cambricon_config.device_indices)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        with safe.block(message="Failed to collect Cambricon MLU metrics", level="debug"):
            mlu_ids = [str(i) for i in self._indices]
            results.extend(self._collect_utilization(mlu_ids))
            results.extend(self._collect_memory(mlu_ids))
            results.extend(self._collect_temperature(mlu_ids))
            results.extend(self._collect_power(mlu_ids))

        return results

    def _collect_utilization(self, mlu_ids: List[str]) -> List[CollectResult]:
        results: List[CollectResult] = []
        for mid in mlu_ids:
            results.append((f"mlu.{mid}.pct", math.nan))
        with safe.block(message="Failed to collect MLU utilization", level="debug"):
            output = subprocess.run(["cnmon", "info", "-u"], capture_output=True, text=True).stdout
            index = 0
            for line in output.split("\n"):
                if "mlu average" in line.lower() and index < len(mlu_ids):
                    util_str = line.split(":")[-1].replace("%", "").strip()
                    try:
                        results[index] = (f"mlu.{mlu_ids[index]}.pct", float(util_str))
                    except ValueError:
                        pass
                    index += 1
        return results

    def _collect_memory(self, mlu_ids: List[str]) -> List[CollectResult]:
        results: List[CollectResult] = []
        for mid in mlu_ids:
            results.append((f"mlu.{mid}.mem.pct", math.nan))
            results.append((f"mlu.{mid}.mem.value", math.nan))
        with safe.block(message="Failed to collect MLU memory", level="debug"):
            output = subprocess.run(["cnmon", "info", "-m"], capture_output=True, text=True).stdout
            lines = output.split("\n")
            index = 0
            for line_idx, line in enumerate(lines):
                if "physical memory usage" in line.lower() and index < len(mlu_ids):
                    if line_idx + 2 < len(lines) and "used" in lines[line_idx + 2].lower():
                        used_line = lines[line_idx + 2]
                        used_str = used_line.split(":")[-1].replace("MiB", "").strip()
                        try:
                            used_val = float(used_str)
                            mlu_id = mlu_ids[index]
                            total_gb = float(self._mlu_info.get(mlu_id, {}).get("memory", 0))
                            if total_gb > 0:
                                pct_idx = index * 2
                                val_idx = index * 2 + 1
                                results[pct_idx] = (f"mlu.{mlu_id}.mem.pct", used_val / (total_gb * 1024) * 100)
                                results[val_idx] = (f"mlu.{mlu_id}.mem.value", used_val)
                        except ValueError:
                            pass
                        index += 1
        return results

    def _collect_temperature(self, mlu_ids: List[str]) -> List[CollectResult]:
        results: List[CollectResult] = []
        for mid in mlu_ids:
            results.append((f"mlu.{mid}.temp", math.nan))
        with safe.block(message="Failed to collect MLU temperature", level="debug"):
            output = subprocess.run(["cnmon", "info", "-e"], capture_output=True, text=True).stdout
            lines = output.split("\n")
            index = 0
            for line_idx, line in enumerate(lines):
                if "temperature" in line.lower() and index < len(mlu_ids):
                    if line_idx + 2 < len(lines) and "chip" in lines[line_idx + 2].lower():
                        temp_line = lines[line_idx + 2]
                        temp_str = temp_line.split(":")[-1].replace("C", "").strip()
                        try:
                            results[index] = (f"mlu.{mlu_ids[index]}.temp", float(temp_str))
                        except ValueError:
                            pass
                        index += 1
        return results

    def _collect_power(self, mlu_ids: List[str]) -> List[CollectResult]:
        results: List[CollectResult] = []
        for mid in mlu_ids:
            results.append((f"mlu.{mid}.power", math.nan))
        with safe.block(message="Failed to collect MLU power", level="debug"):
            output = subprocess.run(["cnmon", "info", "-p"], capture_output=True, text=True).stdout
            lines = output.split("\n")
            index = 0
            for line_idx, line in enumerate(lines):
                if "power" in line.lower() and index < len(mlu_ids):
                    if line_idx + 1 < len(lines) and "usage" in lines[line_idx + 1].lower():
                        power_line = lines[line_idx + 1]
                        power_str = power_line.split(":")[-1].replace("W", "").strip()
                        try:
                            results[index] = (f"mlu.{mlu_ids[index]}.power", float(power_str))
                        except ValueError:
                            pass
                        index += 1
        return results

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get Cambricon MLU info")
    def get() -> Optional[AcceleratorSnapshot]:
        if platform.system() != "Linux":
            console.debug(
                f"Cambricon MLU detection skipped: cnmon is only supported on Linux, current platform is {platform.system()}"
            )
            return None

        driver, mlu_map = CambriconMLU._map_mlu_raw()
        if not mlu_map:
            console.debug("No Cambricon MLU detected: cnmon info -m returned no MLU devices")
            return None

        devices = []
        for idx, mlu_id in enumerate(sorted(mlu_map.keys(), key=int)):
            info = mlu_map[mlu_id]
            devices.append(
                DeviceSnapshot(
                    index=int(mlu_id),
                    name=f"MLU {mlu_id} {info.get('name', '')}",
                    memory=int(info.get("memory", 0)),
                    memory_unit="GB",
                )
            )

        names = [d.name for d in devices]
        console.debug(f"Detected {len(devices)} Cambricon MLU device(s): {', '.join(names)}")
        return AcceleratorSnapshot(
            vendor="cambricon",
            version=driver,
            devices=devices,
        )

    @staticmethod
    def _map_mlu_raw() -> Tuple[Optional[str], dict]:
        output = ""
        with safe.block(message="Failed to run cnmon info for MLU mapping", level="debug"):
            output = subprocess.run(["cnmon", "info", "-m"], capture_output=True, check=True, text=True).stdout
        mlu_map: Dict[str, dict] = {}
        driver = None
        lines = output.split("\n")[1:]

        mlu_id = None
        for line in lines:
            parts = line.split()
            with safe.block(message="Failed to parse MLU info line", level="debug"):
                if parts[0] == "Card":
                    mlu_id = parts[-1]
                    if mlu_id not in mlu_map:
                        mlu_map[mlu_id] = {}
                if mlu_id is None:
                    continue
                if parts[0] == "Product":
                    mlu_map[mlu_id]["name"] = parts[-1]
                if parts[0] == "Driver" and len(mlu_map[mlu_id]) == 1 and driver is None:
                    driver = parts[-1]
                if parts[0] == "Total" and len(mlu_map[mlu_id]) == 1:
                    mlu_map[mlu_id]["memory"] = str(int(parts[-2]) // 1024)

        return driver, mlu_map
