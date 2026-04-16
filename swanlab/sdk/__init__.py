"""
@author: cunyue
@file: __init__.py
@time: 2026/3/5 14:35
@description: SwanLab SDK，负责SwanLab库的核心指标上传功能
"""

from .cmd.init import init
from .cmd.login import login, login_cli
from .cmd.logout import logout_cli
from .cmd.merge_settings import merge_settings
from .cmd.run import async_log, define_scalar, finish, log, log_audio, log_image, log_text, log_video
from .cmd.verify import verify_cli
from .internal.pkg import console, fs, helper, safe
from .internal.run import Run, clear_run, get_run, has_run, set_run
from .internal.run.config import config
from .internal.run.transforms import Audio, Image, Text, Video
from .internal.settings import Settings
from .protocol import Callback

__all__ = [
    # cmd
    "init",
    "finish",
    "login",
    "log",
    "log_text",
    "log_image",
    "log_audio",
    "log_video",
    "async_log",
    "verify_cli",
    "logout_cli",
    "define_scalar",
    "merge_settings",
    "login_cli",
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
