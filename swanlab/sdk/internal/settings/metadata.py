"""
@author: cunyue
@file: metadata.py
@time: 2026/3/5 17:37
@description: SwanLab 元数据代理、监控、采集配置，涉及：
1. 硬件监控
2. 系统信息采集
3. 终端日志采集
"""

import os
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, DirectoryPath, Field


def get_default_system_drive() -> Path:
    # Windows 系统下，SystemDrive 环境变量默认为 C:，需要手动添加末尾的斜杠
    if os.name == "nt":
        return Path(os.environ.get("SystemDrive", "C:")).resolve()
    return Path("/")


class EnvironmentSettings(BaseModel):
    """启动时一次性采集的系统快照。"""

    hardware: bool = True
    """Controls the collection, reporting, and persistence of static hardware metadata (e.g., GPU model, CPU cores, total memory) at startup.

    Note on interaction with `monitor`:
    If `hardware` is set to False while `monitor.enable` is True, the underlying system hardware information will still be accessed by the monitoring module to compute dynamic metrics (like utilization percentages). 
    However, the static hardware snapshot itself will be explicitly discarded — it will neither be included in the telemetry payload nor saved to local persistent storage.
    """

    runtime: bool = True
    """Controls the logging of the software execution environment. 
    When True, captures details such as the operating system, Python version, hostname, current working directory, and the exact command used to launch the script.
    """

    requirements: bool = True
    """Enables the snapshot of Python dependencies. 
    When True, records the installed pip packages and their exact versions (similar to `pip freeze`) to ensure environment reproducibility.
    """

    conda: bool = False
    """Enables the extraction of Conda environment configurations. 
    When True, records the active Conda environment details and exported dependencies. It defaults to False as fetching Conda metadata may introduce slight overhead during startup.
    """

    git: bool = True
    """Controls the tracking of Git repository metadata. 
    When True, captures the current branch, latest commit hash, and remote URL. This helps tightly link the experiment run to a specific version of your codebase.
    """

    # TODO: There are some gpu/npu specific environment variables that can be collected, such as CUDA_VISIBLE_DEVICES, ROC_VISIBLE_DEVICES, etc.

    model_config = ConfigDict(frozen=True)


class MonitorSettings(BaseModel):
    """周期性硬件指标监控。"""

    enable: bool = True
    """Whether to enable periodic hardware monitoring (CPU usage, GPU utilization, memory, etc.)."""
    interval: int = Field(default=10, ge=5)
    """Hardware monitoring collection interval (seconds)."""
    disk_io_dir: DirectoryPath = Field(default_factory=get_default_system_drive)
    """Directory path for disk I/O monitoring."""

    model_config = ConfigDict(frozen=True)


class ConsoleSettings(BaseModel):
    proxy_type: Literal["all", "stdout", "stderr", "none"] = "all"
    """Terminal log proxy strategy."""
    max_log_length: int = Field(default=1024, ge=500, le=4096)
    """Maximum character length per line for terminal log collection."""

    model_config = ConfigDict(frozen=True)
