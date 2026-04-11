"""
@author: cunyue
@file: __init__.py
@time: 2026/3/5 14:35
@description: SwanLab SDK，负责SwanLab库的核心指标上传功能
"""

from swanlab.sdk.cmd.init import init
from swanlab.sdk.cmd.login import login
from swanlab.sdk.cmd.merge_settings import Settings, merge_settings
from swanlab.sdk.cmd.run import async_log, define_scalar, finish, log, log_audio, log_image, log_text, log_video
from swanlab.sdk.internal.pkg.safe import safe, safe_block
from swanlab.sdk.internal.run import Run, clear_run, get_run, has_run, set_run
from swanlab.sdk.internal.run.config import config

__all__ = [
    "merge_settings",
    "Settings",
    "init",
    "finish",
    "login",
    "log",
    "log_text",
    "log_image",
    "log_audio",
    "log_video",
    "async_log",
    "define_scalar",
    "Run",
    "has_run",
    "get_run",
    "set_run",
    "clear_run",
    "config",
    "safe_block",
    "safe",
]
