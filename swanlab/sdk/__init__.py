"""
@author: cunyue
@file: __init__.py
@time: 2026/3/5 14:35
@description: SwanLab SDK，负责SwanLab库的核心指标上传功能
"""

from swanlab.sdk.internal.run.components.config import config

from .cmd import utils as cmd_utils
from .cmd.init import init
from .cmd.login import login, login_cli
from .cmd.merge_callbacks import merge_callbacks
from .cmd.merge_settings import merge_settings
from .cmd.run import (
    async_log,
    define_scalar,
    finish,
    log,
    log_audio,
    log_echarts,
    log_image,
    log_molecule,
    log_object3d,
    log_text,
    log_video,
)
from .internal import pkg
from .internal.run import Run, clear_run, get_run, has_run, set_run
from .internal.run.transforms import Audio, ECharts, Image, Molecule, Object3D, Text, Video, echarts, plot
from .internal.settings import Settings, settings
from .protocol import Callback

__all__ = [
    # cmd
    "init",
    "finish",
    "login",
    "login_cli",
    "log",
    "log_text",
    "log_image",
    "log_audio",
    "log_video",
    "log_echarts",
    "log_object3d",
    "log_molecule",
    "async_log",
    "define_scalar",
    "merge_settings",
    "merge_callbacks",
    "cmd_utils",
    # run
    "has_run",
    "get_run",
    "set_run",
    "clear_run",
    "config",
    # utils
    "pkg",
    "settings",
    # data
    "Audio",
    "ECharts",
    "Image",
    "Molecule",
    "Object3D",
    "Text",
    "Video",
    "echarts",
    "plot",
    # protocol
    "Callback",
    "Settings",
    "Run",
]
