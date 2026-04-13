"""
@author: cunyue
@file: __init__.py
@time: 2026/3/5 14:35
@description: SwanLab SDK，负责SwanLab库的核心指标上传功能
"""

from .cmd import (
    async_log,
    define_scalar,
    finish,
    init,
    log,
    log_audio,
    log_image,
    log_text,
    log_video,
    login,
    login_interactive,
    logout,
    merge_settings,
    verify,
)
from .internal.pkg import console, fs, helper, safe
from .internal.protocol import Callback
from .internal.run import Run, clear_run, get_run, has_run, set_run
from .internal.run.config import config
from .internal.run.transforms import Audio, Image, Text, Video
from .internal.settings import Settings

__all__ = [
    # cmd
    "init",
    "finish",
    "login",
    "logout",
    "verify",
    "log",
    "log_text",
    "log_image",
    "log_audio",
    "log_video",
    "async_log",
    "define_scalar",
    "merge_settings",
    "login_interactive",
    # run
    "has_run",
    "get_run",
    "set_run",
    "clear_run",
    "config",
    # utils
    "safe",
    "helper",
    "console",
    "fs",
    # data
    "Audio",
    "Image",
    "Text",
    "Video",
    # protocol
    "Callback",
    "Settings",
    "Run",
]
