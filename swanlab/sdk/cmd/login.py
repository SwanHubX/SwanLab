"""
@author: cunyue
@file: login.py
@time: 2026/3/6 22:24
@description: swanlab.login 方法，登录到 SwanLab 平台
"""

import getpass
import sys
from pathlib import Path
from typing import Optional

from rich.text import Text

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.cmd import utils
from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import console, fs, helper, nrc, safe, scope
from swanlab.sdk.internal.settings import ROOT_FOLDER, Settings
from swanlab.sdk.internal.settings import settings as global_settings
from swanlab.sdk.typings.cmd import LoginType
from swanlab.sdk.typings.pkg.client.bootstrap import LoginResponse

__all__ = ["login", "login_cli", "login_raw"]


@with_cmd_lock
@without_run("login")
def login(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: LoginType = False,
    timeout: int = 10,
) -> bool:
    """Authenticate with SwanLab Cloud.

    This function authenticates your environment with SwanLab. If already logged in
    and `relogin` is False, this function does nothing. Call this before `swanlab.init()`
    to use cloud features.

    :param api_key: Your SwanLab API key. If not provided, will attempt to read from
        environment or prompt for input.

    :param relogin: If True, forces re-authentication and overwrites existing credentials.
        Defaults to False.

    :param host: Custom API host URL. If not provided, uses the default SwanLab cloud host.

    :param save: Whether to save the API key locally for future sessions. Defaults to False.

    :param timeout: Network request timeout in seconds. Defaults to 10.

    :return: True if login was successful, False otherwise.

    :raises RuntimeError: If called while a run is active.

    :raises AuthenticationError: If login fails due to invalid credentials or network issues.

    Examples:

        Login with an API key:

        >>> import swanlab
        >>> swanlab.login(api_key="your_api_key_here")
        >>> swanlab.init(mode="cloud")

        Interactive login (prompts for API key):

        >>> import swanlab
        >>> swanlab.login()
        >>> swanlab.init(mode="cloud")

        Force re-login and save credentials:

        >>> import swanlab
        >>> swanlab.login(api_key="new_api_key", relogin=True, save=True)
        >>> swanlab.init(mode="cloud")
    """
    return login_raw(
        api_key=api_key,
        relogin=relogin,
        host=host,
        save=save,
        timeout=timeout,
    )


def login_raw(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: LoginType = False,
    timeout: int = 10,
    wellcome_on_success: bool = True,
    animation: bool = True,
) -> bool:
    # 1. 判断是否允许重新登录
    # 如果已经登录且不需要重新登录，则直接返回
    # 仅当运行时 client 已存在时才视为已登录；本地凭证仅表示可复用，不代表本次会话已完成认证
    already_logged_in = client.exists()
    if already_logged_in and not relogin:
        console.info(
            "You are already logged in. Use",
            Text("`swanlab.login(relogin=True)`", style="bold"),
            "to force relogin.",
            sep=" ",
        )
        return True
    if client.exists():
        client.reset()
    # 2. 获取当前配置
    host = nrc.fmt(host) if host is not None else None
    # 先用入参，入参没有才考虑复用 settings 里的值
    if api_key is None:
        # host 变了，且 .netrc 中存有旧凭证 —— 旧 key 与新 host 不匹配，不能复用
        if host is not None and host != global_settings.api_host and global_settings.api_key is not None:
            raise ValueError(
                f"Stored API key is for '{global_settings.api_host}', but you are logging in to '{host}'. "
                "Please provide an API key for the new host."
            )
        else:
            if global_settings.api_key is None:
                raise ValueError("No API key provided and no stored API key found. Please provide an API key.")
            api_key = global_settings.api_key
    api_host = host or global_settings.api_host
    login_settings = Settings.model_validate({"api_key": api_key, "api_host": api_host, "web_host": host})
    # 3. 进入登录流程
    login_settings.merge_settings({"api_key": api_key})
    with scope.Scope() as s:
        f = utils.with_loading_animation("Waiting for response...")(create_client) if animation else create_client
        f(api_key=api_key, api_host=api_host, timeout=timeout)
        login_resp: Optional[LoginResponse] = s.get("login_resp", None)
        if wellcome_on_success:
            wellcome(login_resp)
        if save:
            nrc_path = nrc.path(Path.cwd() / ROOT_FOLDER) if save == "local" else nrc.path(global_settings.root)
            nrc.write(nrc_path, api_host=api_host, web_host=login_settings.web_host, api_key=api_key)
        # 4. 将登录设置合并到全局配置中
        global_settings.merge_settings(login_settings)
        return True


def create_client(api_key: str, api_host: str, timeout: int = 10):
    return client.new(api_key, api_host, timeout=timeout)


def login_cli(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: LoginType = True,
    timeout: int = 10,
) -> bool:
    """
    带循环输入容错的交互式登录接口。
    主要为 CLI 环境或需要极高容错的终端调用设计。
    当捕获到 AuthenticationError 时，如果环境允许交互，则会无限循环提示用户重新输入 API Key。
    """
    assert save is not False, "login_cli must save credentials locally to support CLI usage"
    # CLI 每次是新进程，需要检查本地凭证判断是否已登录
    nrc_path: Path
    nrc_path = nrc.path(Path.cwd() / ROOT_FOLDER) if save == "local" else nrc.path(global_settings.root)
    already_logged_in = nrc.read(nrc_path) is not None
    if already_logged_in and not relogin:
        console.info(
            "You are already logged in. Use",
            Text("`swanlab login --relogin`", style="bold"),
            "to force relogin.",
            sep=" ",
        )
        return True
    if host is not None:
        host = nrc.fmt(host)
        tmp = Settings(api_host=host, web_host=host)
        api_host = tmp.api_host
        web_host = tmp.web_host
    else:
        api_host = Settings.model_fields["api_host"].default
        web_host = Settings.model_fields["web_host"].default
    count = 0
    base_url = api_host + "/api"
    interactive = global_settings.interactive
    while True:
        if not api_key:
            api_key = prompt_api_key(web_host=web_host, interactive=interactive, again=count > 0)
        try:
            with scope.Scope() as s:
                client.new(api_key, base_url, timeout=timeout)
                login_resp: Optional[LoginResponse] = s.get("login_resp", None)
            wellcome(login_resp)
            write_gitignore = save == "local" and not nrc_path.exists()
            if write_gitignore:
                if not nrc_path.parent.exists():
                    fs.safe_mkdirs(nrc_path.parent)
                utils.append_gitignore(nrc_path.parent)
            nrc.write(nrc_path, api_host=api_host, web_host=web_host, api_key=api_key)
            return True
        except AuthenticationError as e:
            # 如果全局配置禁用了交互模式，直接抛出异常
            if not interactive:
                raise e
            console.error(str(e))
            api_key = None
        except (KeyboardInterrupt, EOFError):
            console.info("\nLogin cancelled by user.")
            return False
        count = count + 1


def prompt_api_key(
    web_host: str,
    interactive: bool,
    tip: str = "Paste an API key from your profile and hit enter, or press 'CTRL + C' to quit",
    again: bool = False,
) -> str:
    """
    让用户在终端安全地输入 API Key
    输入时内容将被隐藏。完整保留了原本的交互文案与 Windows 专属提示

    :param web_host: 当前 Web 主机地址
    :param interactive: 全局配置是否为交互模式
    :param tip: 提示信息
    :param again: 是否为重试模式，重试模式下会略微调整提示信息以区分首次输入与重试输入

    :raises RuntimeError: 如果当前环境不支持交互式输入
    :return: 用户输入的 API Key
    """
    if not interactive:
        raise RuntimeError(
            "API Key not provided and interactive mode is disabled",
            "use `swanlab.login(interactive=True)` or SWANLAB_INTERACTIVE=1 to enable interactive mode.",
        )
    if not helper.is_interactive():
        raise RuntimeError("Cannot prompt for API Key in no-tty environment")
    # 1. 打印获取 API Key 的指引（非重试模式下）
    if not again:
        # 动态拼接当前环境的设置页 URL
        setting_url = f"{web_host}/space/~/settings#development"
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
    with safe.block(message="Failed to read API Key from terminal"):
        try:
            # 隐藏输入内容
            key = getpass.getpass("")
            return key.strip()
        except (KeyboardInterrupt, EOFError):
            # 优雅处理用户按下 Ctrl+C 或 Ctrl+D 退出的情况，替代旧版的 sys.excepthook
            console.print("\n")  # 换行，防止终端提示符错位
            sys.exit(0)


def wellcome(login_resp: Optional[LoginResponse]):
    """
    登录成功后打印欢迎信息
    :param login_resp: 登录响应对象，包含用户信息等数据
    :return:
    """
    assert login_resp is not None, "Login response is missing"
    username = login_resp.get("userInfo", {}).get("username", "unknown")
    console.info("Login successfully. Hi", Text(f"{username}!", "bold"), sep=" ")
