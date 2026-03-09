"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 13:20
@description: SwanLab SDK 控制台日志模块，负责打印日志
当 SWANLAB_DEBUG=true 时，终端输出会同步镜像到诊断日志文件（通过 log 模块）
"""

from rich.console import Console
from rich.markup import escape
from rich.text import Text

from swanlab.sdk.internal.pkg import log as _log
from swanlab.sdk.pkg import helper

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


def _to_plain_text(*args) -> str:
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
    if not _can_log or not helper.env.DEBUG:
        return
    _print_formatted("debug", *args, **kwargs)
    _log.debug("[CONSOLE] %s", _to_plain_text(*args))


def info(*args, **kwargs):
    """发送常规通知"""
    if not _can_log:
        return
    _print_formatted("info", *args, **kwargs)
    _log.info("[CONSOLE] %s", _to_plain_text(*args))


def warning(*args, **kwargs):
    """发生警告"""
    if not _can_log:
        return
    _print_formatted("warning", *args, **kwargs)
    _log.warning("[CONSOLE] %s", _to_plain_text(*args))


def error(*args, **kwargs):
    """发生错误"""
    if not _can_log:
        return
    _print_formatted("error", *args, **kwargs)
    _log.error("[CONSOLE] %s", _to_plain_text(*args))


def critical(*args, **kwargs):
    """致命错误"""
    if not _can_log:
        return
    _print_formatted("critical", *args, **kwargs)

    _log.critical("[CONSOLE] %s", _to_plain_text(*args))
