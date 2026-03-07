"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:53
@description: SwanLab API Key 管理，主要与本地存储相关
"""

import netrc
from typing import Optional

from swanlab.sdk.internal.pkg.settings import get_current_settings

from .helper import get_nrc_path, remove_host_suffix

__all__ = ["get", "save", "exists"]


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

    nrc = netrc.netrc(nrc_path)
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
    """
    current_settings = get_current_settings()
    if host is None:
        host = remove_host_suffix(current_settings.api_url, "/api")
    nrc_path = get_nrc_path()
    if not nrc_path.exists():
        current_settings.save_dir.mkdir(parents=True, exist_ok=True)
        nrc_path.touch()
    nrc = netrc.netrc(nrc_path)
    new_info = (username, "", api_key)
    # 避免重复的写
    info = nrc.authenticators(host)
    if info is None or (info[0], info[2]) != (new_info[0], new_info[2]):
        # 同时只允许存在一个host： https://github.com/SwanHubX/SwanLab/issues/797
        nrc.hosts = {host: new_info}
        with open(nrc_path, "w") as f:
            f.write(repr(nrc))
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
