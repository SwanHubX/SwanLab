"""
@author: cunyue
@file: console.py
@time: 2026/5/10 21:44
@description: 控制台输出设置
"""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class TerminalSettings(BaseModel):
    proxy_type: Literal["all", "stdout", "stderr", "none"] = "all"
    """Terminal log proxy strategy."""
    max_log_length: int = Field(default=1024, ge=500, le=4096)
    """Maximum character length per line for terminal log collection."""

    model_config = ConfigDict(frozen=True)
