"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 13:20
@description: SwanLab SDK 控制台日志模块，负责打印日志
当 SWANLAB_DEBUG=true 时，终端输出会同步镜像到诊断日志文件（通过 log 模块）
"""

import inspect
import os
from datetime import datetime
from typing import Any, Literal, Optional

from rich.console import Console
from rich.text import Text

from swanlab.sdk.internal.pkg import log
from swanlab.sdk.utils import helper

__all__ = ["debug", "info", "warning", "error", "trace", "disable_console", "enable_console"]

_console = Console()
_can_log = True
_name = "swanlab"
_this_file = os.path.abspath(__file__)


def disable_console():
    """
    Disable logging.
    """
    global _can_log
    _can_log = False


def enable_console():
    """
    Enable logging.
    """
    global _can_log
    _can_log = True


def _now() -> str:
    """返回毫秒级时间戳，格式：2026-03-14 10:23:45.124"""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]


def _caller_location() -> str:
    """
    遍历调用栈，返回第一个不属于本模块的帧的 module:func:line。
    用于 loguru 风格的调用方定位。
    """
    for fi in inspect.stack():
        if os.path.abspath(fi.filename) != _this_file:
            module = fi.frame.f_globals.get("__name__", "swanlab.internal")
            return f"{module}:{fi.function}:{fi.lineno}"
    return "swanlab.internal"


def _to_plain_text(*args: Any) -> str:
    """
    将 console 的 *args 转为纯文本字符串，用于写入诊断日志文件。
    去除 rich markup，保留原始语义。
    """
    parts = []
    for arg in args:
        if isinstance(arg, Text):
            parts.append(arg.plain)
        else:
            parts.append(str(arg))
    return " ".join(parts)


def _loguru_build(level_name: str, style: str, *args, location: Optional[str] = None) -> tuple["Text", str]:
    """
    构建 loguru 风格的 Rich Text（终端用）和纯文本字符串（日志文件用），
    两者共享同一份时间戳和调用方位置，确保一致性。
    location 未传时自动从调用栈推导。
    """
    now = _now()
    loc: str = location if location is not None else _caller_location()
    level_padded = level_name.upper().ljust(8)

    line = Text()
    line.append(now, style=style)
    line.append(" | ", style=f"{style} dim")
    line.append(level_padded, style=f"{style} bold")
    line.append(" | ", style=f"{style} dim")
    line.append(loc, style=style)
    line.append(" - ", style=f"{style} dim")
    for i, arg in enumerate(args):
        if i > 0:
            line.append(" ")
        if isinstance(arg, Text):
            line.append_text(arg)
        else:
            line.append(str(arg))

    plain = f"{now} | {level_padded} | {loc} - {_to_plain_text(*args)}"
    return line, plain


def _loguru_print(level_name: str, style: str, *args, **kwargs) -> str:
    """
    构建并打印 loguru 风格日志行，同时返回对应的纯文本字符串供写入日志文件。
    2026-03-14 10:23:45.124 | WARNING  | module:func:line - message
    """
    line, plain = _loguru_build(level_name, style, *args)
    _console.print(line, **kwargs)
    return plain


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
    """发送调试消息（仅 SWANLAB_DEBUG=true 时输出）
    格式：2026-03-14 10:23:45.124 | DEBUG    | module:func:line - message
    """
    if not _can_log or not helper.env.DEBUG:
        return
    plain = _loguru_print("debug", "grey54", *args, **kwargs)
    log.debug(plain)


def info(*args, color: str = "blue", **kwargs):
    """发送常规通知
    格式：swanlab: message
    :param color: 前缀 'swanlab' 的颜色，默认蓝色
    """
    if not _can_log:
        return
    prefix = Text(_name, style=f"{color} bold", no_wrap=True) + Text(":", style="default")
    safe_args = [Text(str(a)) if not isinstance(a, Text) else a for a in args]
    _console.print(prefix, *safe_args, **kwargs)


def warning(*args, **kwargs):
    """发生警告
    格式：2026-03-14 10:23:45.124 | WARNING  | module:func:line - message
    """
    if not _can_log:
        return
    plain = _loguru_print("warning", "yellow", *args, **kwargs)
    log.warning(plain)


def error(*args, **kwargs):
    """发生错误
    格式：2026-03-14 10:23:45.124 | ERROR    | module:func:line - message
    """
    if not _can_log:
        return
    plain = _loguru_print("error", "red", *args, **kwargs)
    log.error(plain)


_LOG_FUNCS = {
    "debug": log.debug,
    "warning": log.warning,
    "error": log.error,
}


def trace(*args, level: Literal["debug", "warning", "error"] = "debug", id: Optional[str] = None):  # noqa: A002
    """仅写入日志文件，不打印到终端。
    格式：2026-03-14 10:23:45.124 | DEBUG    | <id 或调用栈> - message

    :param level: 日志级别，可选 'debug' / 'warning' / 'error'，默认 'debug'
    :param id:    位置标识符；传入时直接使用，省略时自动从调用栈推导
    """
    if not _can_log:
        return
    log_func = _LOG_FUNCS.get(level, log.debug)
    _, plain = _loguru_build(level, "", *args, location=id)
    log_func(plain)
