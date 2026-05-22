"""
@author: cunyue
@file: amd.py
@time: 2026/3/31 01:57
@description: AMD ROCm GPU 信息采集模块

检测原理：
- Linux 优先读取 amdgpu 驱动暴露的 sysfs 文件，避免 rocm-smi 在部分环境中不稳定。
- 设备扫描格式：遍历 `/sys/class/drm/card*`，读取 `device/vendor`，AMD vendor id 为 `0x1002`；
  总显存读取 `device/mem_info_vram_total`（bytes）。驱动版本读取 `/sys/module/amdgpu/version`。
- Linux 动态指标格式：`gpu_busy_percent` 为 0-100；`mem_info_vram_used/total` 为 bytes；
  `hwmon/hwmon*/temp1_input` 为 milli-degree Celsius；`power1_average` 或 `power1_input` 为 microwatt。
- Windows 兜底使用 ROCm `hipinfo` 文本输出，按 `device# <idx>` 分块，解析 `Name:`、
  `totalGlobalMem: <num> GB`、`memInfo.total/free: <num> GB`。
- 如果 sysfs 扫描不到设备，会尝试 `rocm-smi --showproductname --json` 的 `cardX -> Card Series` JSON。
"""

import glob
import json
import math
import os
import platform
import re
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


class AMDGPU(AcceleratorProtocol):
    def __init__(self, shim: SystemShim):
        super().__init__(shim)
        amd_config = next(a for a in shim.accelerators if a.vendor == "rocm")
        self._indices = amd_config.device_indices
        self._system = platform.system()
        self._hwmon_paths: Dict[int, Optional[str]] = {}
        if self._system == "Linux":
            for idx in self._indices:
                base = f"/sys/class/drm/card{idx}/device"
                hwmons = glob.glob(os.path.join(base, "hwmon", "hwmon*"))
                self._hwmon_paths[idx] = hwmons[0] if hwmons else None

    @classmethod
    @safe.decorator(level="debug", message="Failed to initialize AMD monitor")
    def new(cls, shim: SystemShim) -> Optional[Tuple["AMDGPU", SystemScalars]]:
        amd_config = next((a for a in shim.accelerators if a.vendor == "rocm"), None)
        if amd_config is None:
            console.debug("No AMD monitor config in shim")
            return None
        self = cls(shim)
        max_mem_mb = 0
        for gpu_id in amd_config.device_indices:
            if self._system == "Linux":
                with safe.block(message=f"Failed to read AMD GPU {gpu_id} memory info", level="debug"):
                    base = f"/sys/class/drm/card{gpu_id}/device"
                    with open(os.path.join(base, "mem_info_vram_total"), "r") as f:
                        total_bytes = int(f.read().strip())
                    total_mb = total_bytes // (1024 * 1024)
                    if total_mb > max_mem_mb:
                        max_mem_mb = total_mb

        scalars: SystemScalars = []
        for color_idx, idx in enumerate(amd_config.device_indices):
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

            mem_pct = SystemScalar(
                key=f"gpu.{idx}.mem.pct",
                name=f"GPU {idx}",
                chart_name="GPU Memory Allocated (%)",
                y_min=0,
                y_max=100,
                color=color,
            )
            scalars.append(mem_pct)

            mem_value = SystemScalar(
                key=f"gpu.{idx}.mem.value",
                name=f"GPU {idx}",
                chart_name="GPU Memory Allocated (MB)",
                y_min=0,
                y_max=max_mem_mb,
                color=color,
            )
            scalars.append(mem_value)

            temp = SystemScalar(
                key=f"gpu.{idx}.temp",
                name=f"GPU {idx}",
                chart_name="GPU Temperature (°C)",
                color=color,
            )
            scalars.append(temp)

            power = SystemScalar(
                key=f"gpu.{idx}.power",
                name=f"GPU {idx}",
                chart_name="GPU Power (W)",
                y_min=0,
                color=color,
            )
            scalars.append(power)

        console.debug(f"AMD monitor initialized with {len(amd_config.device_indices)} device(s)")
        return self, scalars

    def collect(self) -> List[CollectResult]:
        with safe.block(message="Failed to collect AMD GPU metrics", level="debug"):
            if self._system == "Linux":
                return self._collect_linux_sysfs()
            else:
                return self._collect_windows()
        return []

    def _collect_linux_sysfs(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        for idx in self._indices:
            base = f"/sys/class/drm/card{idx}/device"

            results.append((f"gpu.{idx}.pct", math.nan))
            with safe.block(message=f"Failed to collect AMD GPU {idx} utilization", level="debug"):
                with open(os.path.join(base, "gpu_busy_percent"), "r") as f:
                    util = float(f.read().strip())
                results[-1] = (f"gpu.{idx}.pct", util)

            mem_pct = math.nan
            mem_used_mb = math.nan
            with safe.block(message=f"Failed to collect AMD GPU {idx} memory", level="debug"):
                with open(os.path.join(base, "mem_info_vram_used"), "r") as f:
                    used_bytes = int(f.read().strip())
                with open(os.path.join(base, "mem_info_vram_total"), "r") as f:
                    total_bytes = int(f.read().strip())
                if total_bytes > 0:
                    mem_pct = (used_bytes / total_bytes) * 100
                    mem_used_mb = used_bytes / (1024 * 1024)
            results.append((f"gpu.{idx}.mem.pct", mem_pct))
            results.append((f"gpu.{idx}.mem.value", mem_used_mb))

            temp_val = math.nan
            power_val = math.nan
            with safe.block(message=f"Failed to collect AMD GPU {idx} temperature/power", level="debug"):
                hwmon_path = self._hwmon_paths.get(idx)
                if hwmon_path is not None:
                    if os.path.exists(os.path.join(hwmon_path, "temp1_input")):
                        with open(os.path.join(hwmon_path, "temp1_input"), "r") as f:
                            temp_val = float(f.read().strip()) / 1000.0
                    p_file = os.path.join(hwmon_path, "power1_average")
                    if not os.path.exists(p_file):
                        p_file = os.path.join(hwmon_path, "power1_input")
                    if os.path.exists(p_file):
                        with open(p_file, "r") as f:
                            power_val = float(f.read().strip()) / 1000000.0
            results.append((f"gpu.{idx}.temp", temp_val))
            results.append((f"gpu.{idx}.power", power_val))

        return results

    def _collect_windows(self) -> List[CollectResult]:
        results: List[CollectResult] = []
        output = ""
        with safe.block(message="Failed to run hipinfo for AMD GPU collection", level="debug"):
            output = subprocess.run(["hipinfo"], capture_output=True, text=True).stdout

        devices_content: Dict[str, str] = {}
        parts = re.split(r"device#\s+(\d+)", output)
        if len(parts) >= 2:
            for i in range(1, len(parts), 2):
                devices_content[parts[i].strip()] = parts[i + 1]
        elif output:
            devices_content["0"] = output

        for idx in self._indices:
            content = devices_content.get(str(idx), "")

            mem_total = self._extract_float(content, r"memInfo\.total:\s*([\d\.]+)\s*GB")
            mem_free = self._extract_float(content, r"memInfo\.free:\s*([\d\.]+)\s*GB")

            mem_pct = math.nan
            mem_used_mb = math.nan
            if mem_total is not None and mem_free is not None and mem_total > 0:
                mem_used = mem_total - mem_free
                mem_used_mb = mem_used * 1024
                mem_pct = (mem_used / mem_total) * 100

            results.append((f"gpu.{idx}.pct", math.nan))
            results.append((f"gpu.{idx}.mem.pct", mem_pct))
            results.append((f"gpu.{idx}.mem.value", mem_used_mb))
            results.append((f"gpu.{idx}.temp", math.nan))
            results.append((f"gpu.{idx}.power", math.nan))

        return results

    @staticmethod
    @safe.decorator(level="debug", message="Failed to extract float from hipinfo output")
    def _extract_float(content: str, regex: str) -> Optional[float]:
        match = re.search(regex, content, re.IGNORECASE)
        if match:
            return float(match.group(1))
        return None

    @staticmethod
    @safe.decorator(level="debug", message="Failed to get AMD GPU info")
    def get() -> Optional[AcceleratorSnapshot]:
        system_name = platform.system()
        if system_name not in ["Linux", "Windows"]:
            console.debug(
                f"AMD GPU detection skipped: only Linux/Windows are supported, current platform is {system_name}"
            )
            return None

        if system_name == "Linux":
            driver, gpu_map = AMDGPU._map_amd_gpu_linux()
        else:
            driver, gpu_map = AMDGPU._map_amd_gpu_windows()

        if not gpu_map:
            console.debug("No AMD GPU detected: sysfs/hipinfo/rocm-smi did not return any AMD devices")
            return None

        devices = []
        for gpu_id, info in gpu_map.items():
            mem_str = info.get("memory", "0GB")
            digits = re.findall(r"([\d\.]+)", mem_str)
            mem_val = int(float(digits[0])) if digits else 0
            unit = "GB" if "GB" in mem_str else "MB"
            devices.append(
                DeviceSnapshot(index=int(gpu_id), name=info.get("name", "AMD GPU"), memory=mem_val, memory_unit=unit)
            )

        names = [d.name for d in devices]
        console.debug(f"Detected {len(devices)} AMD GPU(s): {', '.join(names)}")
        return AcceleratorSnapshot(vendor="rocm", version=driver, devices=devices)

    @staticmethod
    def _map_amd_gpu_linux() -> Tuple[Optional[str], dict]:
        driver_version = "Unknown"
        with safe.block(message="Failed to get AMD GPU driver version", level="debug"):
            if os.path.exists("/sys/module/amdgpu/version"):
                with open("/sys/module/amdgpu/version", "r") as f:
                    driver_version = f.read().strip()

        gpu_map: dict = {}
        with safe.block(message="Failed to map AMD GPUs via sysfs", level="debug"):
            cards = glob.glob("/sys/class/drm/card*")
            amd_cards = []
            for card_path in cards:
                if not re.search(r"card\d+$", card_path):
                    continue
                vendor_path = os.path.join(card_path, "device/vendor")
                if os.path.exists(vendor_path):
                    with open(vendor_path, "r") as f:
                        vendor = f.read().strip()
                    if "0x1002" in vendor:
                        idx = card_path.split("card")[-1]
                        amd_cards.append(idx)
            amd_cards.sort(key=lambda x: int(x))

            for idx in amd_cards:
                mem_str = "0GB"
                mem_total_path = f"/sys/class/drm/card{idx}/device/mem_info_vram_total"
                if os.path.exists(mem_total_path):
                    with open(mem_total_path, "r") as f:
                        mem_bytes = int(f.read().strip())
                        mem_gb = round(mem_bytes / (1024**3))
                        mem_str = f"{mem_gb}GB"
                gpu_map[idx] = {"name": "AMD GPU", "memory": mem_str}

        if not gpu_map:
            with safe.block(message="Failed to map AMD GPUs via rocm-smi fallback", level="debug"):
                data = json.loads(
                    subprocess.run(["rocm-smi", "--showproductname", "--json"], capture_output=True, text=True).stdout
                )
                for k, v in data.items():
                    idx = k.replace("card", "")
                    gpu_map[idx] = {"name": v.get("Card Series", "AMD GPU"), "memory": "0GB"}
        return driver_version, gpu_map

    @staticmethod
    def _map_amd_gpu_windows() -> Tuple[Optional[str], dict]:
        driver_version = "Unknown"
        gpu_map: dict = {}
        with safe.block(message="Failed to map AMD GPUs via hipinfo on Windows", level="debug"):
            output = subprocess.run(["hipinfo"], capture_output=True, text=True, check=True).stdout
            devices = re.split(r"device#\s+(\d+)", output)
            if len(devices) < 2:
                gpu_map.update(AMDGPU._parse_hipinfo_text_block("0", output))
            else:
                for i in range(1, len(devices), 2):
                    dev_id = devices[i].strip()
                    gpu_map.update(AMDGPU._parse_hipinfo_text_block(dev_id, devices[i + 1]))
        return driver_version, gpu_map

    @staticmethod
    def _parse_hipinfo_text_block(dev_id: str, content: str) -> dict:
        name = "AMD GPU"
        memory = "0GB"
        name_match = re.search(r"Name:\s*(.+)", content)
        if name_match:
            name = name_match.group(1).strip()
        mem_match = re.search(r"totalGlobalMem:\s*([\d\.]+)\s*GB", content)
        if mem_match:
            memory = f"{mem_match.group(1)}GB"
        return {dev_id: {"name": name, "memory": memory}}
