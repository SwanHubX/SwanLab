"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 01:18
@description: cmd 模块工具函数
"""

import sys
from functools import wraps
from pathlib import Path

from swanlab.sdk.internal.pkg import console, nrc
from swanlab.sdk.internal.settings import settings
from swanlab.sdk.typings.cmd import LoginType

__all__ = ["get_nrc_path", "with_loading_animation", "prompt_masked"]


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


def prompt_masked(mask: str = "*") -> str:
    """在终端逐字符读取一行输入，屏幕上以遮罩符号回显。

    基于 ``pwinput`` 实现跨平台的星号遮罩输入。在 pwinput 基础上补齐两项行为：
    Ctrl+C / Ctrl+D 抛出与 ``getpass`` 一致的异常（pwinput 原生会吞掉），
    以及禁用括号化粘贴模式以防粘贴标记污染输入。

    :param mask: 遮罩符号，默认 "*"
    :return: 用户输入的原始字符串（未 strip）
    :raises KeyboardInterrupt: 用户按下 Ctrl+C
    :raises EOFError: 用户按下 Ctrl+D / Ctrl+Z
    """
    import pwinput

    _orig_getch = pwinput.getch

    def _getch() -> str:
        ch = _orig_getch()
        o = ord(ch)
        if o == 3:  # Ctrl+C
            raise KeyboardInterrupt
        if o in (4, 26):  # Ctrl+D (POSIX) / Ctrl+Z (Windows)
            raise EOFError
        return ch

    # 禁用括号化粘贴模式，防止 \x1b[200~...\x1b[201~ 标记被当作可打印字符注入输入
    sys.stdout.write("\x1b[?2004l")
    sys.stdout.flush()
    pwinput.getch = _getch
    try:
        return pwinput.pwinput(prompt="", mask=mask)
    finally:
        pwinput.getch = _orig_getch
        sys.stdout.write("\x1b[?2004h")
        sys.stdout.flush()
