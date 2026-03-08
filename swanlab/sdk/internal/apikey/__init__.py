"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:53
@description: SwanLab API Key 管理，会根据当前settings获取、保存API Key
"""

from typing import Optional

from swanlab.sdk.internal.pkg.netrc import get_nrc_path, remove_host_suffix, write_netrc
from swanlab.sdk.internal.pkg.settings import get_current_settings

__all__ = ["get", "save", "exists"]


def save(username: str, api_key: str, host: Optional[str] = None):
    """
    保存API Key到本地存储，并更新当前settings中的API Key
    """
    current_settings = get_current_settings()

    if host is None:
        host = remove_host_suffix(current_settings.api_host, "/api")

    nrc_path = get_nrc_path(current_settings.root)

    # 调用底层工具写入凭证
    write_netrc(nrc_path=nrc_path, host=host, username=username, password=api_key)

    # 同步更新运行时 Settings
    current_settings.merge_settings({"api_key": api_key})


def exists() -> bool:
    """
    检查当前的API Key是否存在
    """
    return get_current_settings().api_key is not None


def get() -> str:
    """
    获取当前的API Key，如果不存在则报错
    """
    api_key = get_current_settings().api_key
    if api_key is None:
        raise FileNotFoundError("The API key file or target host does not exist. Please check your API Key.")
    return api_key


def prompt() -> str:
    """
    提示用户输入API Key，如果交互模式未启用，则抛出异常
    输入API Key时会使用password保护，有部分终端不支持，此时依然抛出异常
    """
    # current_settings = get_current_settings()
    # if not current_settings.interactive:
    #     raise RuntimeError(
    #         "API Key not provided and interactive mode is disabled",
    #         "use `swanlab.login(interactive=True)` or SWANLAB_INTERACTIVE=1 to enable interactive mode.",
    #     )
    ...
