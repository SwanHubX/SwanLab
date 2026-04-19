"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 13:20
@description: SwanLab SDK 控制台日志模块，负责打印日志
当 SWANLAB_DEBUG=true 时，终端输出会同步镜像到诊断日志文件（通过 log 模块）
"""

import inspect
import os
import sys
from datetime import datetime
from typing import Any, Literal, Optional

from rich.console import Console
from rich.text import Text

from .. import helper
from . import log

__all__ = ["debug", "info", "warning", "error", "trace", "c", "init", "reset"]

c = Console()
_name = "swanlab"
_this_file = os.path.abspath(__file__)


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


def _loguru_build(
    level_name: str,
    style: str,
    *args,
    location: Optional[str] = None,
    console_args=None,
    file_args=None,
) -> tuple["Text", str]:
    """
    构建 loguru 风格的 Rich Text（终端用）和纯文本字符串（日志文件用），
    两者共享同一份时间戳和调用方位置，确保一致性。
    location 未传时自动从调用栈推导。

    - args: 默认同时用于终端和文件
    - console_args: 若传入，则仅用于终端消息体
    - file_args: 若传入，则仅用于文件消息体
    """
    now = _now()
    loc: str = location if location is not None else _caller_location()
    level_padded = level_name.upper().ljust(8)

    if console_args is None:
        console_args = args
    if file_args is None:
        file_args = args

    line = Text()
    line.append(now, style=style)
    line.append(" | ", style=f"{style} dim")
    line.append(level_padded, style=f"{style} bold")
    line.append(" | ", style=f"{style} dim")
    line.append(loc, style=style)
    line.append(" - ", style=f"{style} dim")

    for i, arg in enumerate(console_args):
        if i > 0:
            line.append(" ")
        if isinstance(arg, Text):
            line.append_text(arg)
        else:
            line.append(str(arg))

    plain = f"{now} | {level_padded} | {loc} - {_to_plain_text(*file_args)}"
    return line, plain


def _loguru_print(level_name: str, style: str, *args, **kwargs) -> str:
    """
    构建并打印 loguru 风格日志行，同时返回对应的纯文本字符串供写入日志文件。
    2026-03-14 10:23:45.124 | DEBUG  | module:func:line - message
    """
    console_args = kwargs.pop("console_args", None)
    file_args = kwargs.pop("file_args", None)
    location = kwargs.pop("location", None)

    line, plain = _loguru_build(
        level_name,
        style,
        *args,
        location=location,
        console_args=console_args,
        file_args=file_args,
    )
    c.print(line, **kwargs)
    return plain


# -----------------------------------------------------------------------------
# 导出的模块级日志函数
# -----------------------------------------------------------------------------


# noinspection PyShadowingBuiltins
def print(*args, **kwargs):  # noqa: A001
    """发送普通消息"""
    c.print(*args, **kwargs)


def debug(*args, write_to_file: bool = True, **kwargs):
    """发送调试消息（仅 SWANLAB_DEBUG=true 时输出）
    格式：2026-03-14 10:23:45.124 | DEBUG    | module:func:line - message
    """
    if not helper.DEBUG:
        return
    plain = _loguru_print("debug", "grey54", *args, **kwargs)
    if write_to_file:
        log.debug(plain)


def info(*args, color: str = "blue", **kwargs):
    """发送常规通知
    格式：swanlab: message
    :param color: 前缀 'swanlab' 的颜色，默认蓝色
    """
    prefix = Text(_name, style=f"{color} bold", no_wrap=True) + Text(":", style="default")
    safe_args = [Text(str(a)) if not isinstance(a, Text) else a for a in args]
    c.print(prefix, *safe_args, **kwargs)
    if helper.DEBUG:
        _, plain = _loguru_build("info", "", *args)
        log.info(plain)


def warning(*args, **kwargs):
    """发生警告
    格式：2026-03-14 10:23:45.124 | WARNING  | module:func:line - message
    """
    plain = _loguru_print("warning", "yellow", *args, **kwargs)
    log.warning(plain)


def error(*args, write_to_file: bool = True, **kwargs):
    """发生错误
    格式：2026-03-14 10:23:45.124 | ERROR    | module:func:line - message
    """
    plain = _loguru_print("error", "red", *args, **kwargs)
    if write_to_file:
        log.error(plain)


def trace(
    *args,
    max_frames: int = 2,
    write_to_file: bool = True,
    level_name: Literal["error", "debug"] = "error",
    **kwargs,
):
    """
    打印当前异常栈，是 error/debug 的变体：在消息末尾追加异常栈信息。
    必须在 except 块内调用，否则无操作。

    - 终端：显示最后 max_frames 层栈帧 + 异常类型 + 异常信息
    - 日志文件：显示完整错误栈

    Args:
        *args: 前缀消息
        max_frames: 终端显示的最大栈帧数，默认 2
        write_to_file: 是否写入日志文件，默认 True
        level_name: 日志级别，"error" 或 "debug"，默认 "error"
    """
    import traceback

    exc_type, exc_value, exc_tb = sys.exc_info()
    if exc_type is None:
        return None

    # 获取完整的错误栈（用于日志文件）
    full_tb = "".join(traceback.format_exception(exc_type, exc_value, exc_tb)).rstrip()
    # 获取最后 N 层栈帧（用于终端）
    frames = traceback.extract_tb(exc_tb)
    if not frames:
        short_tb = f"{exc_type.__name__}: {exc_value}"
    else:
        selected = frames[-max_frames:]
        stack = " -> ".join(f"{frame.name}:{frame.lineno}" for frame in selected)
        short_tb = f"{stack} -> {exc_type.__name__}: {exc_value}"

    prefix = _to_plain_text(*args).strip()
    console_msg = f"{prefix}: {short_tb}" if prefix else short_tb
    file_msg = f"{prefix}: {full_tb}" if prefix else full_tb

    log_fn = error if level_name == "error" else debug
    log_fn(console_args=(console_msg,), file_args=(file_msg,), write_to_file=write_to_file, **kwargs)


# -----------------------------------------------------------------------------
# 日志文件管理（委托给 log 模块）
# -----------------------------------------------------------------------------


def init(bind_to=None) -> None:
    """初始化日志模块。bind_to 为日志目录时绑定文件，为 None 时禁用持久化。仅应被 Run 生命周期调用。"""
    log.init(bind_to)


def reset() -> None:
    """重置日志模块到初始状态。仅应被 Run 生命周期和测试 teardown 调用。"""
    log.reset()
