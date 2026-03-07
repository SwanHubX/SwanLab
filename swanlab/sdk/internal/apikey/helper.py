"""
@author: cunyue
@file: helper.py
@time: 2026/3/7 14:56
@description: API Key管理工具函数
"""

from pathlib import Path

from swanlab.sdk.internal.pkg.settings import get_current_settings


def remove_host_suffix(host: str, *suffixes: str) -> str:
    """
    移除host的后缀
    :param host: 待处理的host
    :param suffixes: 要移除的后缀列表
    :return: 处理后的host
    """
    host = host.rstrip()
    if not suffixes:
        return host
    for suffix in suffixes:
        if not suffix:
            continue
        if host.endswith(suffix):
            return host[: -len(suffix)]
    return host


def get_nrc_path() -> Path:
    """
    获取netrc文件路径，并不保证文件、目录存在
    """
    current_settings = get_current_settings()
    return current_settings.root / ".netrc"
