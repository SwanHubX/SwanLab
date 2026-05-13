"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 01:18
@description: cmd 模块工具函数
"""

from functools import wraps
from pathlib import Path

from swanlab.sdk.internal.pkg import console, nrc
from swanlab.sdk.internal.settings import settings
from swanlab.sdk.typings.cmd import LoginType

__all__ = ["get_nrc_path", "with_loading_animation"]


def get_nrc_path(save: LoginType) -> Path:
    """根据登录类型获取 NRC 文件路径"""
    assert save, "LoginType cannot be False when getting NRC path"
    return nrc.path(settings.get_pwd_config_dir() if save == "local" else settings.get_user_config_dir())


def with_loading_animation(message: str = "Initializing SwanLab...", spinner_name: str = "dots"):
    """
    一个使用 rich 包装的加载动画装饰器
    :param message: 动画旁边显示的文字
    :param spinner_name: 动画样式（如 'dots', 'bouncingBar', 'point' 等）
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # 使用 rich 的 status 作为上下文管理器包裹函数的执行
            with console.c.status(message, spinner=spinner_name):
                return func(*args, **kwargs)

        return wrapper

    return decorator
