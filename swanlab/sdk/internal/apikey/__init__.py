"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:53
@description: SwanLab API Key 管理，主要与本地存储相关
"""

import netrc
import os
import stat
from pathlib import Path
from typing import Optional

from swanlab.sdk.internal.pkg.settings import get_current_settings

from .helper import get_nrc_path, remove_host_suffix

__all__ = ["get", "save", "exists"]


def create_nrc(path: Path) -> netrc.netrc:
    try:
        return netrc.netrc(path)
    except (netrc.NetrcParseError, IsADirectoryError):
        raise IOError(
            f"Failed to access or parse netrc file at {path}. "
            f"Please check if the path is a directory, has incorrect permissions, "
            f"or contains syntax errors."
        )


def _get_or_none() -> Optional[str]:
    """
    底层核心获取逻辑，安全地尝试获取 API Key，失败返回 None
    """
    current_settings = get_current_settings()
    if current_settings.api_key:
        return current_settings.api_key

    nrc_path = get_nrc_path()
    if not nrc_path.exists():
        return None

    try:
        nrc = netrc.netrc(nrc_path)
    except (netrc.NetrcParseError, IsADirectoryError):
        raise IOError(
            f"Failed to access or parse netrc file at {nrc_path}. "
            f"Please check if the path is a directory, has incorrect permissions, "
            f"or contains syntax errors."
        )
    host = remove_host_suffix(current_settings.api_url, "/api")
    info = nrc.authenticators(host)

    # 向下兼容逻辑
    if info is None:
        info = nrc.authenticators(host + "/api")
        if info is not None:
            nrc.hosts = {host: info}
            with open(nrc_path, "w") as f:
                f.write(repr(nrc))

    if info is None:
        return None

    return info[2]


def save(username: str, api_key: str, host: Optional[str] = None):
    """
    保存API Key到本地存储，并更新当前settings中的API Key
    :param username: 保存的用户名，默认为user，可选择存储为前端网页ip或者域名
    :param api_key: 保存的API Key
    :param host: 保存的host
    :raise IOError: 如果保存过程中出现错误
    :raise OSError: 如果文件权限设置失败或者无法写入文件
    """
    current_settings = get_current_settings()
    if host is None:
        host = remove_host_suffix(current_settings.api_url, "/api")
    nrc_path = get_nrc_path()

    # 1. 核心防御：处理路径冲突（如果是目录则无法继续）
    if nrc_path.exists() and nrc_path.is_dir():
        raise IOError(
            f"The path {nrc_path} is a directory, but a file is expected. "
            f"Please remove the directory or change the SWANLAB_ROOT settings."
        )

    # 2. 确保父目录存在
    nrc_path.parent.mkdir(parents=True, exist_ok=True)

    # 3. 文件初始化与权限控制
    if not nrc_path.exists():
        try:
            # 创建文件并限制权限
            nrc_path.touch(mode=0o600)
        except OSError as e:
            raise OSError(f"Failed to create netrc file at {nrc_path}: {e}. Check your permissions.")
    else:
        os.chmod(nrc_path, stat.S_IRUSR | stat.S_IWUSR)

    # 4. 尝试加载或重置 netrc 对象
    try:
        # 即使文件为空，部分版本的 netrc() 也会报错，所以加一层 try
        nrc = netrc.netrc(nrc_path)
    except (netrc.NetrcParseError, IOError, Exception):
        nrc_path.write_text("")
        nrc = netrc.netrc(nrc_path)

    # 5. 写入逻辑
    new_info = (username, "", api_key)
    info = nrc.authenticators(host)

    if info is None or (info[0], info[2]) != (new_info[0], new_info[2]):
        nrc.hosts = {host: new_info}
        try:
            with open(nrc_path, "w") as f:
                f.write(repr(nrc))
            os.chmod(nrc_path, stat.S_IRUSR | stat.S_IWUSR)
        except OSError as e:
            raise OSError(f"Failed to write API Key to {nrc_path}: {e}")

    current_settings.merge_settings({"api_key": api_key})


def exists() -> bool:
    """
    检查当前的API Key是否存在
    """
    return _get_or_none() is not None


def get() -> str:
    """
    获取当前的API Key，如果不存在则报错
    """
    api_key = _get_or_none()
    if api_key is None:
        # 统一抛出错误，你可以考虑以后换成自定义的 APIKeyNotFoundError
        raise FileNotFoundError("The API key file or target host does not exist. Please check your API Key.")
    return api_key


def prompt() -> str:
    """
    提示用户输入API Key
    """
    ...
