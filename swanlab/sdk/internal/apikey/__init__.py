"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:53
@description: SwanLab API Key 管理，会根据当前settings获取、保存API Key
"""

import getpass
import sys
from typing import Optional

from rich.text import Text

from swanlab.sdk.internal.pkg.netrc import get_nrc_path, remove_host_suffix, write_netrc
from swanlab.sdk.internal.pkg.settings import get_current_settings
from swanlab.sdk.pkg import console, helper

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


def prompt(
    tip: str = "Paste an API key from your profile and hit enter, or press 'CTRL + C' to quit",
    again: bool = False,
) -> str:
    """
    让用户在终端安全地输入 API Key。
    输入时内容将被隐藏。完整保留了原本的交互文案与 Windows 专属提示。

    :param tip: 提示信息
    :param again: 是否是重新输入，如果是，则不显示获取 Key 的链接 URL
    :raises RuntimeError: 如果当前环境不支持交互式输入
    :return: 用户输入的 API Key
    """
    current_settings = get_current_settings()
    if not current_settings.interactive:
        raise RuntimeError(
            "API Key not provided and interactive mode is disabled",
            "use `swanlab.login(interactive=True)` or SWANLAB_INTERACTIVE=1 to enable interactive mode.",
        )
    if not helper.env.is_interactive():
        raise RuntimeError("Cannot prompt for API Key in no-tty environment")
    web_host = current_settings.web_host
    # 1. 打印获取 API Key 的指引（非重试模式下）
    if not again:
        # 动态拼接当前环境的设置页 URL
        setting_url = f"{web_host.rstrip('/')}/space/~/settings#development"
        console.info("You can find your API key at:", Text(setting_url, style="yellow"))

    # 2. 拼接输入提示语
    prompt_text = tip

    # 针对 Windows 环境的专属粘贴提示
    if sys.platform == "win32":
        prompt_text += (
            "\nOn Windows, use [yellow]Ctrl + Shift + V[/yellow] or [yellow]right-click[/yellow] to paste the API key"
        )

    prompt_text += ": "

    # 先使用 console 打印提示，因为 getpass() 原生不支持 Rich 的颜色标签渲染
    console.print(prompt_text, end="")

    # 强制刷新输出缓冲区，确保提示语立刻显示
    sys.stdout.flush()

    # 3. 安全读取用户输入
    try:
        # 隐藏输入内容
        key = getpass.getpass("")
        return key.strip()
    except (KeyboardInterrupt, EOFError):
        # 优雅处理用户按下 Ctrl+C 或 Ctrl+D 退出的情况，替代旧版的 sys.excepthook
        console.print("\n")  # 换行，防止终端提示符错位
        sys.exit(0)
    except Exception as e:
        raise RuntimeError(f"Failed to read API Key from terminal: {e}")
