"""
@author: cunyue
@file: kunlunxin.py
@time: 2026/3/31 01:58
@description: 昆仑芯 Kunlunxin XPU 信息采集模块

检测原理：
- 仅 Linux 支持，依赖昆仑芯 `xpu-smi` 文本输出。
- 驱动版本使用 `xpu-smi -q`，解析包含 `Driver Version` 的行，取冒号后的值。
- 设备和动态指标均使用 `xpu-smi -m` 的空白分隔表格。当前解析依赖字段位置：
  第 2 列为 XPU ID（index 1），第 5 列为温度（index 4），第 9 列为功耗（index 8），
  第 18 列为已用显存 MB（index 17），第 19 列为总显存 MB（index 18），第 20 列为利用率（index 19），
  第 22/23 列拼接为设备名称（index 21/22）。
- 显存使用率按 `used_mb / total_mb * 100` 计算；旧版未提供独立显存使用量曲线，当前沿用该行为。
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


class KunlunxinXPU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        kunlun_config = next(a for a in shim.accelerators if a.vendor == "kunlunxin")
        self._indices = kunlun_config.device_indices
        driver, xpu_map = KunlunxinXPU._map_xpu()
        self._xpu_map = xpu_map

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize Kunlunxin XPU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["KunlunxinXPU", SystemScalars]]:
        kunlun_config = next((a for a in shim.accelerators if a.vendor == "kunlunxin"), None)
        if kunlun_config is None:
            console.debug("No Kunlunxin XPU monitor config in shim")
            return None
        self = cls(shim)

        scalars: SystemScalars = []
        for color_idx, idx in enumerate(self._indices):
            color = generate_color(color_idx)

            scalars.append(
                SystemScalar(
                    key=f"xpu.{idx}.pct",
                    name=f"XPU {idx}",
                    chart_name="XPU Utilization (%)",
                    y_min=0,
                    y_max=100,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"xpu.{idx}.mem.pct",
                    name=f"XPU {idx}",
                    chart_name="XPU Memory Allocated (%)",
                    y_min=0,
                    y_max=100,
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"xpu.{idx}.temp",
                    name=f"XPU {idx}",
                    chart_name="XPU Temperature (°C)",
                    color=color,
                )
            )
            scalars.append(
                SystemScalar(
                    key=f"xpu.{idx}.power",
                    name=f"XPU {idx}",
                    chart_name="XPU Power (W)",
                    y_min=0,
                    color=color,
                )
            )

        console.debug(f"Kunlunxin XPU monitor initialized with {len(self._indices)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        with safe.block(message="Failed to collect Kunlunxin XPU metrics", level="debug"):
            metrics_by_id: Dict[int, Dict[str, float]] = {idx: {} for idx in self._indices}
            output = subprocess.run(["xpu-smi", "-m"], capture_output=True, text=True).stdout
            for line in output.split("\n"):
                parts = line.split()
                if len(parts) < 23:
                    continue
                xpu_id_str = parts[1]
                if not xpu_id_str.isdigit() or int(xpu_id_str) not in self._indices:
                    continue
                xpu_id = int(xpu_id_str)
                with safe.block(message="Failed to parse Kunlunxin XPU row data", level="debug"):
                    util_str = parts[19]
                    mem_used_str = parts[17]
                    temp_str = parts[4]
                    power_str = parts[8]
                    try:
                        metrics_by_id[xpu_id]["pct"] = float(util_str)
                    except ValueError:
                        pass
                    xpu_info = self._xpu_map.get(xpu_id_str, {})
                    total_mb = int(xpu_info.get("memory", 0)) * 1024
                    try:
                        if total_mb > 0:
                            metrics_by_id[xpu_id]["mem.pct"] = float(mem_used_str) / total_mb * 100
                    except ValueError:
                        pass
                    try:
                        metrics_by_id[xpu_id]["temp"] = float(temp_str)
                    except ValueError:
                        pass
                    try:
                        metrics_by_id[xpu_id]["power"] = float(power_str)
                    except ValueError:
                        pass

            results: List[CollectResult] = []
            for idx in self._indices:
                m = metrics_by_id[idx]
                results.append((f"xpu.{idx}.pct", m.get("pct", math.nan)))
                results.append((f"xpu.{idx}.mem.pct", m.get("mem.pct", math.nan)))
                results.append((f"xpu.{idx}.temp", m.get("temp", math.nan)))
                results.append((f"xpu.{idx}.power", m.get("power", math.nan)))
            return results
        return []

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get Kunlunxin XPU info")
    def get() -> Optional[AcceleratorSnapshot]:
        if platform.system() != "Linux":
            console.debug(
                f"Kunlunxin XPU detection skipped: xpu-smi is only supported on Linux, current platform is {platform.system()}"
            )
            return None

        driver, xpu_map = KunlunxinXPU._map_xpu()
        if not xpu_map:
            console.debug("No Kunlunxin XPU detected: xpu-smi -m returned no parseable XPU rows")
            return None

        devices = []
        for xpu_id in sorted(xpu_map.keys(), key=int):
            info = xpu_map[xpu_id]
            devices.append(
                DeviceSnapshot(
                    index=int(xpu_id),
                    name=info.get("name", "Kunlunxin XPU"),
                    memory=int(info.get("memory", 0)),
                    memory_unit="GB",
                )
            )

        names = [d.name for d in devices]
        console.debug(f"Detected {len(devices)} Kunlunxin XPU device(s): {', '.join(names)}")
        return AcceleratorSnapshot(
            vendor="kunlunxin",
            version=driver,
            devices=devices,
        )

    @staticmethod
    def _map_xpu() -> Tuple[Optional[str], dict]:
        output = subprocess.run(["xpu-smi", "-m"], capture_output=True, check=True, text=True).stdout

        driver = KunlunxinXPU._get_xpu_driver()
        xpu_map: Dict[str, dict] = {}

        for line in output.split("\n"):
            parts = line.split()
            with safe.block(message="Failed to parse Kunlunxin XPU info row", level="debug"):
                if len(parts) >= 23:
                    xpu_id = parts[1]
                    name = f"{parts[21]} {parts[22]}".strip().strip("()")
                    memory = int(parts[18]) // 1024
                    xpu_map[xpu_id] = {"name": name, "memory": str(memory)}

        return driver, xpu_map

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get Kunlunxin XPU driver version")
    def _get_xpu_driver() -> Optional[str]:
        output = subprocess.run(["xpu-smi", "-q"], capture_output=True, check=True, text=True).stdout
        for line in output.split("\n"):
            if "Driver Version" in line:
                return line.split(":")[-1].strip()
        return None
