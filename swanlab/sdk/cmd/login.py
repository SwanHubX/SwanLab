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
from urllib.parse import urlparse

from rich.text import Text

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal.context import RunConfig, RunContext, use_context
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import console, helper, safe, scope
from swanlab.sdk.internal.settings import Settings, settings
from swanlab.sdk.typings.pkg.client.bootstrap import LoginResponse

__all__ = ["login", "login_interactive", "login_raw"]


@with_cmd_lock
@without_run("login")
def login(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: bool = False,
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
    return login_raw(api_key=api_key, relogin=relogin, host=host, save=save, timeout=timeout)


def login_interactive(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: bool = False,
    timeout: int = 10,
) -> bool:
    """
    带循环输入容错的交互式登录接口。
    主要为 CLI 环境或需要极高容错的终端调用设计。
    当捕获到 AuthenticationError 时，如果环境允许交互，则会无限循环提示用户重新输入 API Key。
    """
    # CLI 每次是新进程，client.exists() 必为 False，需要检查本地凭证判断是否已登录
    if settings.api_key is not None and not relogin:
        console.info(
            "You are already logged in. Use",
            Text("`swanlab login --relogin`", style="bold"),
            "to force relogin.",
            sep=" ",
        )
        return True
    try:
        # 首次尝试登录，复用原子接口
        return login_raw(api_key=api_key, relogin=relogin, host=host, save=save, timeout=timeout)
    except AuthenticationError as e:
        # 如果全局配置禁用了交互模式，直接抛出异常
        if not settings.interactive:
            raise e
        console.error(str(e))

    # 进入容错重试循环
    while True:
        try:
            # 重新要求用户输入新的 Key
            new_key = prompt_api_key()
            return login_raw(api_key=new_key, relogin=relogin, host=host, save=save, timeout=timeout)
        except AuthenticationError as e:
            console.error(str(e))
        except (KeyboardInterrupt, EOFError):
            console.info("\nLogin cancelled by user.")
            return False


def login_raw(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: bool = False,
    timeout: int = 10,
) -> bool:
    # 1. 判断是否允许重新登录
    # 如果已经登录且不需要重新登录，则直接返回
    # 仅当运行时 client 已存在时才视为已登录；本地凭证仅表示可复用，不代表本次会话已完成认证
    already_logged_in = client.exists()
    if not relogin and already_logged_in:
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
    # 如果提供了host且没有http前缀添加https://前缀
    if host is not None:
        host = host.strip().rstrip("/")
        if not host.startswith(("http://", "https://")):
            host = f"https://{host}"
        parsed = urlparse(host)
        host = f"{parsed.scheme}://{parsed.netloc}"
    # 先用入参，入参没有才考虑复用 settings 里的值
    if api_key is None:
        # host 变了，且 .netrc 中存有旧凭证 —— 旧 key 与新 host 不匹配，不能复用
        if host is not None and host != settings.api_host and settings.api_key is not None:
            console.warning(
                f"Stored API key is for '{settings.api_host}', but you are logging in to '{host}'. "
                "Please provide an API key for the new host."
            )
            # api_key 保持 None，后续会走 prompt
        else:
            api_key = settings.api_key
    api_host = host or settings.api_host
    # 防止用户通过 login(host="api.swanlab.cn") 强行覆盖 web_host
    if host and "api.swanlab.cn" in host:
        web_host = settings.web_host  # 保持官方默认 web_host
    else:
        web_host = host or settings.web_host
    login_settings = Settings.model_validate({"api_key": api_key, "api_host": api_host, "web_host": web_host})
    # 临时使用运行上下文，在结束后清除
    fake_run_dir = Path.cwd()
    with use_context(RunContext(config=RunConfig(settings=login_settings, run_dir=fake_run_dir))) as ctx:
        # 如果 API Key 不存在，则提示用户输入
        if api_key is None:
            api_key = prompt_api_key(ctx=ctx)
        # 3. 进入登录流程
        ctx.config.settings.merge_settings({"api_key": api_key})
        with scope.Scope() as s:
            create_client(ctx, timeout=timeout)
            assert client.exists(), "Failed to create client"
            login_resp: Optional[LoginResponse] = s.get("login_resp", None)
            if login_resp is None:
                # 认证失败但 client 对象已被 new() 创建，必须清理否则后续重试会报 "client already exists"
                if client.exists():
                    client.reset()
                raise AuthenticationError("Failed to login, please check your API Key or network connection.")
            username = login_resp.get("userInfo", {}).get("username", "unknown")
            console.info("Login successfully. Hi", Text(f"{username}!", "bold"), sep=" ")
            if save:
                raise NotImplementedError("暂未实现")
        # 4. 将登录设置合并到全局配置中
        settings.merge_settings(ctx.config.settings)
        return True


@helper.with_loading_animation()
def create_client(ctx: RunContext, timeout: int = 10):
    assert ctx.config.settings.api_key is not None, "API Key not provided"
    return client.new(ctx.config.settings.api_key, ctx.config.settings.api_host, timeout=timeout)


def prompt_api_key(
    tip: str = "Paste an API key from your profile and hit enter, or press 'CTRL + C' to quit",
    again: bool = False,
    ctx: Optional[RunContext] = None,
) -> str:
    """
    让用户在终端安全地输入 API Key。
    输入时内容将被隐藏。完整保留了原本的交互文案与 Windows 专属提示。

    :param tip: 提示信息
    :param again: 是否是重新输入，如果是，则不显示获取 Key 的链接 URL
    :param ctx: 当前运行上下文，如果提供则使用上下文中的配置信息，否则使用全局配置信息

    :raises RuntimeError: 如果当前环境不支持交互式输入
    :return: 用户输入的 API Key
    """
    current_settings = settings if ctx is None else ctx.config.settings
    if not current_settings.interactive:
        raise RuntimeError(
            "API Key not provided and interactive mode is disabled",
            "use `swanlab.login(interactive=True)` or SWANLAB_INTERACTIVE=1 to enable interactive mode.",
        )
    if not helper.is_interactive():
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
    with safe.block(message="Failed to read API Key from terminal"):
        try:
            # 隐藏输入内容
            key = getpass.getpass("")
            return key.strip()
        except (KeyboardInterrupt, EOFError):
            # 优雅处理用户按下 Ctrl+C 或 Ctrl+D 退出的情况，替代旧版的 sys.excepthook
            console.print("\n")  # 换行，防止终端提示符错位
            sys.exit(0)
