"""
@author: cunyue
@file: __init__.py
@time: 2026/3/5 14:35
@description: SwanLab SDK，负责SwanLab库的核心指标上传功能
"""

from .cmd.init import init
from .cmd.login import login
from .cmd.logout import logout
from .cmd.merge_settings import Settings, merge_settings
from .cmd.run import async_log, define_scalar, finish, log, log_audio, log_image, log_text, log_video
from .cmd.verify import verify
from .internal.pkg.helper import get_swanlab_version
from .internal.pkg.safe import block, decorator
from .internal.protocol import Callback
from .internal.run import Run, clear_run, get_run, has_run, set_run
from .internal.run.config import config
from .internal.run.transforms import Audio, Image, Text, Video

__all__ = [
    "Callback",
    "get_swanlab_version",
    "Audio",
    "Image",
    "Text",
    "Video",
    "merge_settings",
    "Settings",
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
    "Run",
    "has_run",
    "get_run",
    "set_run",
    "clear_run",
    "config",
    "block",
    "decorator",
]
