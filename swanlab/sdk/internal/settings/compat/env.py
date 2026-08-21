"""
遗留环境变量读取工具：统一剥离 shell 脚本注入的成对引号。

某些 shell 脚本（eval、命令替换、/etc/environment）注入的环境变量会保留字面引号字符，
例如 ``SWANLAB_API_KEY='"abc123"'`` 的实际值为 ``"abc123"``（含引号）。
所有遗留环境变量的读取都应通过 :func:`getenv` 完成，避免引号进入 Settings 解析。
"""

import os
from typing import Optional, overload

__all__ = ["getenv", "strip_env_quotes"]


def strip_env_quotes(value: str) -> str:
    """剥离环境变量值首尾成对的同类引号。

    仅当首尾字符相同且均为双引号或单引号时才剥离，
    避免误伤仅单侧带引号或不匹配的情况。
    """
    if not isinstance(value, str) or len(value) < 2:
        return value
    first, last = value[0], value[-1]
    if first == last and first in ('"', "'"):
        return value[1:-1]
    return value


@overload
def getenv(key: str) -> Optional[str]: ...


@overload
def getenv(key: str, default: None) -> Optional[str]: ...


@overload
def getenv(key: str, default: str) -> str: ...


def getenv(key: str, default: Optional[str] = None) -> Optional[str]:
    """读取环境变量并剥离首尾成对引号。"""
    value = os.environ.get(key, default)
    if isinstance(value, str):
        return strip_env_quotes(value)
    return value
