"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 13:20
@description: SwanLab SDK 控制台日志模块，负责打印日志
"""

import os

from rich.console import Console
from rich.markup import escape
from rich.text import Text

# 设计上 console 独立于其他模块，包括settings。但是我们又希望有一个环境变量去控制是否打印debug信息，所以这里额外绑定一个环境变量
DEBUG = os.getenv("SWANLAB_DEBUG", "false").lower() in ["true", "1", "yes", "on"]


__all__ = ["debug", "info", "warning", "error", "critical", "disable_log", "enable_log"]


_console = Console()
_can_log = True
_name = "swanlab"

# 暂时保留权重设计
_config = {
    "debug": (10, Text(_name, style="grey54 bold", no_wrap=True) + Text(":", style="default")),
    "info": (20, Text(_name, style="blue bold", no_wrap=True) + Text(":", style="default")),
    "warning": (30, Text(_name, style="yellow bold", no_wrap=True) + Text(":", style="default")),
    "error": (40, Text(_name, style="red bold", no_wrap=True) + Text(":", style="default")),
    "critical": (50, Text(_name, style="red bold", no_wrap=True) + Text(":", style="default")),
}


def disable_log():
    """
    Disable logging.
    """
    global _can_log
    _can_log = False


def enable_log():
    """
    Enable logging.
    """
    global _can_log
    _can_log = True


def _print_formatted(level_name: str, *args, **kwargs):
    """
    纯粹的打印执行逻辑，剥离了所有的校验，减少判断开销
    """
    _, prefix = _config[level_name]
    prefix_text = prefix.copy()

    if kwargs.get("sep") == "":
        prefix_text.append(" ")
    # 对 args 里的每一个字符串进行转义
    safe_args = [escape(str(arg)) if isinstance(arg, str) else arg for arg in args]
    _console.print(prefix_text, *safe_args, **kwargs)


# -----------------------------------------------------------------------------
# 导出的模块级日志函数
# -----------------------------------------------------------------------------


# noinspection PyShadowingBuiltins
def print(*args, **kwargs):  # noqa: A001
    """发送普通消息"""
    if not _can_log:
        return
    _console.print(*args, **kwargs)


def debug(*args, **kwargs):
    """发送调试消息"""
    if not _can_log or not DEBUG:
        return
    _print_formatted("debug", *args, **kwargs)


def info(*args, **kwargs):
    """发送常规通知"""
    if not _can_log:
        return
    _print_formatted("info", *args, **kwargs)


def warning(*args, **kwargs):
    """发生警告"""
    if not _can_log:
        return
    _print_formatted("warning", *args, **kwargs)


def error(*args, **kwargs):
    """发生错误"""
    if not _can_log:
        return
    _print_formatted("error", *args, **kwargs)


def critical(*args, **kwargs):
    """致命错误"""
    if not _can_log:
        return
    _print_formatted("critical", *args, **kwargs)
