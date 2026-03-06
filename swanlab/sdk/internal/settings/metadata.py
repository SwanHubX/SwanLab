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


class HardwareSettings(BaseModel):
    enable: bool = True
    """Whether to collect machine hardware information."""
    monitor: bool = True
    """Whether to collect machine hardware monitoring information."""
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


class EnvSettings(BaseModel):
    runtime: bool = True
    """Whether to collect machine runtime information."""
    requirements: bool = True
    """Whether to collect Python environment (requirements) information."""
    conda: bool = False
    """Whether to collect Conda environment information."""
    git: bool = True
    """Whether to collect Git information."""

    model_config = ConfigDict(frozen=True)
