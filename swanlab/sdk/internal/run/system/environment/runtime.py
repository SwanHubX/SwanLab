"""
@author: cunyue
@file: runtime.py
@time: 2026/3/30
@description: 运行时信息采集
"""

import os
import platform
import socket
import sys
from typing import Optional


def get_os() -> str:
    """获取操作系统平台"""
    return platform.platform()


def get_os_pretty() -> Optional[str]:
    """获取操作系统友好名称"""
    return platform.freedesktop_os_release().get("PRETTY_NAME")


def get_hostname() -> str:
    """获取主机名"""
    return socket.gethostname()


def get_pid() -> int:
    """获取进程 ID"""
    return os.getpid()


def get_cwd() -> str:
    """获取当前工作目录"""
    return os.getcwd()


def get_python_version() -> str:
    """获取 Python 版本"""
    return platform.python_version()


def get_python_verbose() -> str:
    """获取 Python 详细版本信息"""
    return sys.version


def get_python_executable() -> str:
    """获取 Python 可执行文件路径"""
    return sys.executable


def get_command() -> str:
    """获取当前执行命令"""
    if platform.system() == "Linux":
        cmdline = _get_command_linux()
        if cmdline:
            return cmdline
    return " ".join(sys.argv)


def _get_command_linux() -> str:
    """Linux 下获取当前执行命令"""
    with open("/proc/self/cmdline", "rb") as f:
        content = f.read()
        if not content:
            raise RuntimeError("Cannot read /proc/self/cmdline")

        # 使用 rb 读取并用 null 字符分割，处理编码更安全
        # cmdline 以 \0 分隔，最后一个通常也是 \0
        args = content.decode("utf-8", errors="replace").split("\x00")
        # 过滤掉空的参数并重新组合
        return " ".join([arg for arg in args if arg]).strip()
