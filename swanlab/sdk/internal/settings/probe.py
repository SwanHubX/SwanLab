"""
@author: cunyue
@file: probe.py
@time: 2026/5/10 21:40
@description: 探针模块设置
"""

import os
from pathlib import Path

from pydantic import BaseModel, ConfigDict, DirectoryPath, Field


def get_default_system_drive() -> Path:
    # Windows 系统下，SystemDrive 环境变量默认为 C:，需要手动添加末尾的斜杠
    if os.name == "nt":
        return Path(os.environ.get("SystemDrive", "C:")).resolve()
    return Path("/")


class ProbeSettings(BaseModel):
    hardware: bool = True
    """Controls the collection, reporting, and persistence of static hardware metadata (e.g., GPU model, CPU cores, total memory) at startup.

    Note on interaction with `monitor`:
    If `hardware` is set to False while `monitor` is True, the underlying system hardware information will still be accessed by the monitoring module to compute dynamic metrics (like utilization percentages). 
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

    swanlab: bool = True
    """Controls the tracking of SwanLab metadata. 
    When True, captures the SwanLab version, the SwanLab run directory and so on.
    """

    # TODO: There are some gpu/npu specific environment variables that can be collected, such as CUDA_VISIBLE_DEVICES, ROC_VISIBLE_DEVICES, etc.

    monitor: bool = True
    """Whether to enable periodic hardware monitoring (CPU usage, GPU utilization, memory, etc.).
    """

    monitor_interval: int = Field(default=10, ge=5)
    """Periodic hardware monitoring interval (seconds)."""

    monitor_disk_dir: DirectoryPath = Field(default_factory=get_default_system_drive)
    """Disk I/O monitoring directory."""

    model_config = ConfigDict(frozen=True)
