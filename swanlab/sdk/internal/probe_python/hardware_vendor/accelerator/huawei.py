"""
@author: cunyue
@file: huawei.py
@time: 2026/3/31 01:58
@description: 华为昇腾 Ascend NPU 信息采集模块

检测原理：
- 仅 Linux 支持，先检查 `/dev/davinci*` 设备文件，再调用华为官方 `npu-smi` 工具。
- 静态映射使用 `npu-smi info -m`，文本表格从第二行开始解析：
  `NPU ID  Chip ID  Chip Logic ID  Chip Name...`；其中 Chip Logic ID 非数字的 MCU 行会被过滤。
- 驱动版本使用 `npu-smi -v`，解析冒号后的版本号。
- CANN 版本读取 `/usr/local/Ascend/ascend-toolkit/latest/{arch}-linux/ascend_toolkit_install.info`
  中的 `version=<value>`。
- 动态指标使用 `npu-smi info -t usages/temp/power -i <npu_id> -c <chip_id>`，解析
  `Aicore Usage Rate`、`HBM Usage Rate`、`HBM Capacity`、温度和功耗冒号后的数值。
"""

import math
import os
import platform
import subprocess
from typing import List, Optional, Tuple

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


class AscendNPU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        ascend_config = next(a for a in shim.accelerators if a.vendor == "ascend")
        self._indices = ascend_config.device_indices
        self._npu_map = AscendNPU._map_npu_raw()
        self._chips: List[Tuple[str, str]] = []
        flat_idx = 0
        for npu_id in sorted(self._npu_map.keys(), key=int):
            for chip_id in sorted(self._npu_map[npu_id].keys(), key=int):
                if flat_idx in self._indices:
                    self._chips.append((npu_id, chip_id))
                flat_idx += 1

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize Ascend NPU monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["AscendNPU", SystemScalars]]:
        ascend_config = next((a for a in shim.accelerators if a.vendor == "ascend"), None)
        if ascend_config is None:
            console.debug("No Ascend NPU monitor config in shim")
            return None
        self = cls(shim)

        hbm_value = 0
        for npu_id, chip_id in self._chips:
            usage = AscendNPU._get_chip_usage(npu_id, chip_id)
            if usage:
                hbm = int(usage.get("hbm", 0))
                if hbm > hbm_value:
                    hbm_value = hbm
        max_hbm_mb = hbm_value * 1024

        scalars: SystemScalars = []
        for color_idx, (npu_id, chip_id) in enumerate(self._chips):
            label = f"{npu_id}-{chip_id}"
            color = generate_color(color_idx)

            util = SystemScalar(
                key=f"npu.{label}.ptc",
                name=f"NPU {label}",
                chart_name="NPU Utilization (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(util)

            hbm_rate = SystemScalar(
                key=f"npu.{label}.mem.ptc",
                name=f"NPU {label}",
                chart_name="NPU Memory Allocated (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(hbm_rate)

            hbm_value_scalar = SystemScalar(
                key=f"npu.{label}.mem.value",
                name=f"NPU {label}",
                chart_name="NPU Memory Allocated (MB)",
                y_min=0,
                y_max=max_hbm_mb,
                color=color,
            )
            scalars.append(hbm_value_scalar)

            temp = SystemScalar(
                key=f"npu.{label}.temp",
                name=f"NPU {label}",
                chart_name="NPU Temperature (°C)",
                color=color,
            )
            scalars.append(temp)

            power = SystemScalar(
                key=f"npu.{label}.power",
                name=f"NPU {label}",
                chart_name="NPU Power Usage (W)",
                color=color,
            )
            scalars.append(power)

        console.debug(f"Ascend NPU monitor initialized with {len(self._chips)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        with safe.block(message="Failed to collect Ascend NPU metrics", level="debug"):
            results: List[CollectResult] = []
            for npu_id, chip_id in self._chips:
                label = f"{npu_id}-{chip_id}"
                results.extend(self._collect_usage(npu_id, chip_id, label))
                results.append(self._collect_temp(npu_id, chip_id, label))
                results.append(self._collect_power(npu_id, chip_id, label))
            return results
        return []

    def _collect_usage(self, npu_id: str, chip_id: str, label: str) -> List[CollectResult]:
        util_val = math.nan
        hbm_val = math.nan
        hbm_mb = math.nan
        with safe.block(message="Failed to collect Ascend NPU usage", level="debug"):
            output = subprocess.run(
                ["npu-smi", "info", "-t", "usages", "-i", npu_id, "-c", chip_id],
                capture_output=True,
                text=True,
            ).stdout
            for line in output.split("\n"):
                if "aicore usage rate" in line.lower():
                    util_str = line.split(":")[-1].strip()
                    if util_str.isdigit():
                        util_val = float(util_str)
                if "hbm usage rate" in line.lower():
                    hbm_str = line.split(":")[-1].strip()
                    if hbm_str.isdigit():
                        hbm_val = float(hbm_str)
        return [
            (f"npu.{label}.ptc", util_val),
            (f"npu.{label}.mem.ptc", hbm_val),
            (f"npu.{label}.mem.value", hbm_mb),
        ]

    def _collect_temp(self, npu_id: str, chip_id: str, label: str) -> CollectResult:
        temp_val = math.nan
        with safe.block(message="Failed to collect Ascend NPU temperature", level="debug"):
            output = subprocess.run(
                ["npu-smi", "info", "-t", "temp", "-i", npu_id, "-c", chip_id],
                capture_output=True,
                text=True,
            ).stdout.strip()
            temp_val = float(output.split(":")[-1].strip())
        return (f"npu.{label}.temp", temp_val)

    def _collect_power(self, npu_id: str, chip_id: str, label: str) -> CollectResult:
        power_val = math.nan
        with safe.block(message="Failed to collect Ascend NPU power", level="debug"):
            output = subprocess.run(
                ["npu-smi", "info", "-t", "power", "-i", npu_id, "-c", chip_id],
                capture_output=True,
                text=True,
            ).stdout.strip()
            power_val = float(output.split(":")[-1].strip())
        return (f"npu.{label}.power", power_val)

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get Ascend NPU info")
    def get() -> Optional[AcceleratorSnapshot]:
        if platform.system() != "Linux":
            console.debug(
                f"Ascend NPU detection skipped: npu-smi is only supported on Linux, current platform is {platform.system()}"
            )
            return None
        if not any(f.startswith("davinci") for f in os.listdir("/dev")):
            console.debug("Ascend NPU detection skipped: no /dev/davinci* device files found")
            return None

        driver = AscendNPU._get_version()
        cann = AscendNPU._get_cann_version()
        npu_map = AscendNPU._map_npu_raw()

        devices = []
        flat_idx = 0
        for npu_id in sorted(npu_map.keys(), key=int):
            for chip_id in sorted(npu_map[npu_id].keys(), key=int):
                chip_info = npu_map[npu_id][chip_id]
                usage = AscendNPU._get_chip_usage(npu_id, chip_id)
                hbm = int(usage.get("hbm", 0)) if usage else 0
                devices.append(
                    DeviceSnapshot(
                        index=flat_idx,
                        name=f"NPU {npu_id}-{chip_id} {chip_info.get('name', '')}",
                        memory=hbm,
                        memory_unit="GB",
                    )
                )
                flat_idx += 1

        if not devices:
            console.debug("No Ascend NPU detected: npu-smi info -m returned no compute chips")
            return None

        console.debug(f"Detected {len(devices)} Ascend NPU device(s): {', '.join(d.name for d in devices)}")
        return AcceleratorSnapshot(
            vendor="ascend",
            version=driver,
            cann_version=cann,
            devices=devices,
        )

    @staticmethod
    def _get_version() -> str:
        result = subprocess.run(["npu-smi", "-v"], capture_output=True, check=True, text=True)
        return result.stdout.split(":")[-1].strip()

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get CANN version")
    def _get_cann_version() -> Optional[str]:
        arch = platform.machine()
        if arch not in ("aarch64", "x86_64"):
            return None
        path = f"/usr/local/Ascend/ascend-toolkit/latest/{arch}-linux/ascend_toolkit_install.info"
        cmd = f"grep -m 1 '^version=' {path} | cut -d'=' -f2"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return proc.stdout.strip() if proc.returncode == 0 else None

    @staticmethod
    def _map_npu_raw() -> dict:
        output = subprocess.run(["npu-smi", "info", "-m"], capture_output=True, check=True, text=True).stdout
        npu_map: dict = {}
        for line in output.split("\n")[1:]:
            parts = line.split()
            if len(parts) < 4:
                continue
            npu_id, chip_id, chip_logic_id, *chip_name = parts
            if not chip_logic_id.isdigit():
                continue
            chip_name_str = " ".join(chip_name)
            if npu_id not in npu_map:
                npu_map[npu_id] = {}
            npu_map[npu_id][chip_id] = {"id": chip_logic_id, "name": chip_name_str}
        return npu_map

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get Ascend NPU chip usage")
    def _get_chip_usage(npu_id: str, chip_id: str) -> Optional[dict]:
        output = subprocess.run(
            ["npu-smi", "info", "-t", "usages", "-i", npu_id, "-c", chip_id],
            capture_output=True,
            text=True,
        ).stdout
        for line in output.split("\n"):
            if "hbm capacity" in line.lower():
                hbm = line.split(":")[-1].strip()
                if hbm.isdigit():
                    return {"hbm": str(round(int(hbm) / 1024))}
        return None
