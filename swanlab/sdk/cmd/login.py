"""
@author: cunyue
@file: login.py
@time: 2026/3/6 22:24
@description: swanlab.login 方法，登录到 SwanLab 平台。
凭据解析契约：
1. 显式提供 API host 时，由 API host 推导 Web host；否则从 Settings 读取 API host 和 Web host。
2. 显式提供 API key 时优先使用该值；否则仅在未提供 API host，或 API host 与 Settings 一致时复用
   Settings 中的 API key。API host 发生变化时不得复用旧 key，且无可用 key 时登录失败。
"""

import sys
from typing import Optional, Tuple

from rich.text import Text

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.cmd import utils
from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal.pkg import console, helper, nrc, safe
from swanlab.sdk.internal.pkg.client import verify_api_key
from swanlab.sdk.internal.settings import Settings, create_settings, resolve_hosts, set_global_settings
from swanlab.sdk.typings.cmd import LoginType
from swanlab.sdk.typings.pkg.client.bootstrap import UserProfile

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

    This function authenticates your environment with SwanLab by verifying the
    provided (or stored) API key online. Every call performs an online check so
    revoked keys fail immediately. Call this before `swanlab.init()` to use cloud
    features; the runtime client itself is created and managed by the core during
    `swanlab.init()`.

    :param api_key: Your SwanLab API key. If not provided, will attempt to read from
        environment or prompt for input.
    :param relogin: Kept for signature compatibility only. Login always verifies
        online, so this flag no longer changes behavior. Defaults to False.
    :param host: Custom API host URL. If not provided, uses the default SwanLab cloud host.
    :param save: Whether to save the API key locally for future sessions. Defaults to False.
    :param timeout: Network request timeout in seconds. Defaults to 10.
    :return: True if login was successful, False otherwise.
    :raises RuntimeError: If called while a run is active.
    :raises AuthenticationError: If credentials are rejected or unavailable.
    :raises RequestException: If a network or HTTP request fails.
    :raises RuntimeError: If the backend does not support ApiKey authentication.
    :raises ValueError: If the authentication response is malformed.

    Examples:

        Login with an API key:

        >>> import swanlab
        >>> swanlab.login(api_key="your_api_key_here")
        >>> swanlab.init(mode="online")

        Force re-login and save credentials:

        >>> import swanlab
        >>> swanlab.login(api_key="new_api_key", relogin=True, save=True)
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
    animation: bool = True,
    print_welcome: bool = True,
) -> bool:
    """显式登录：每次调用都在线验证凭证并获取用户档案，随后关闭临时 client。

    - 不驻留任何运行时 client 单例：online run 的 client 由 Core 在启动时基于 Proto 凭据自建。
    - `relogin` 仅保留为接口签名兼容：登录始终在线验证，不再承载"跳过/强制"分支。
    - `save` 独立控制是否将凭证持久化到本地 netrc 文件。
    """
    # 1. 获取当前 api key、web host、api host 的配置
    current_settings = create_settings()
    host = nrc.fmt(host) if host is not None else None
    # 先用入参，入参没有才考虑复用 settings 里的值
    # 这里存在问题，如果用户输入了host但是没有输入api key，并且当前本地存在旧的api key，那么会导致host和api key不匹配，登录失败
    # 所以这里额外增加api key不存在的时候的判断
    if api_key is None:
        # host 变了，且 .netrc 中存有旧凭证 —— 旧 key 与新 host 不匹配，不能复用
        if host is not None and host != current_settings.api_host and current_settings.api_key is not None:
            raise AuthenticationError(
                f"Stored API key is for '{current_settings.api_host}', but you are logging in to '{host}'. "
                "Please provide an API key for the new host."
            )
        else:
            if current_settings.api_key is None:
                raise AuthenticationError("No API key provided and no stored API key found. Please provide an API key.")
            api_key = current_settings.api_key
    api_host, web_host = validate_host(host, current_settings)
    login_settings = Settings.model_validate({"api_key": api_key, "api_host": api_host, "web_host": web_host})
    # 至此，api_key、api_host、web_host 都已经确定，且 login_settings 已经准备好

    # 2. 临时认证：在线验证凭证并取得 profile，随后立即关闭临时 client（无副作用，不创建运行时单例）
    f = utils.with_loading_animation("Waiting for response...")(verify_api_key) if animation else verify_api_key
    profile = f(api_key=api_key, base_url=api_host, timeout=timeout)
    if print_welcome:
        welcome(api_host, profile)
    # 3. 持久化（由 save 单独控制）并将登录设置合并到全局配置
    if save:
        nrc_path = utils.get_nrc_path(save=save)
        nrc.write(nrc_path, api_host=api_host, web_host=login_settings.web_host, api_key=api_key)
    current_settings.merge_settings(login_settings)
    set_global_settings(current_settings)
    return True


def login_cli(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: LoginType = True,
    timeout: int = 10,
) -> bool:
    """
    带循环输入容错的交互式登录接口。
    本地凭证仅用于选择 API key，每次调用都会在线校验。relogin 保留为签名兼容。
    当捕获到 AuthenticationError 时，如果环境允许交互，则会无限循环提示用户重新输入 API Key。
    """
    assert save is not False, "login_cli must save credentials locally to support CLI usage"
    nrc_path = utils.get_nrc_path(save)
    current_settings = create_settings()
    api_host, web_host = validate_host(host, current_settings)
    stored = nrc.read(nrc_path)
    if api_key is None and stored is not None and stored[1] == api_host:
        api_key = stored[0]

    count = 0
    interactive = current_settings.interactive
    while True:
        if not api_key:
            api_key = prompt_api_key(web_host=web_host, interactive=interactive, again=count > 0)
        try:
            # 临时认证：验证凭证并取得 profile，随后关闭临时 client，不驻留运行时单例
            profile = verify_api_key(api_key=api_key, base_url=api_host, timeout=timeout)
            welcome(api_host, profile)
            # 如果存储当前目录，添加gitignore文件
            # mkdir_and_append_gitignore 自动判断是否为空文件夹，如果是则写入
            if save == "local":
                helper.mkdir_and_append_gitignore(nrc_path.parent)
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


def validate_host(host: Optional[str], settings: Settings) -> Tuple[str, str]:
    """
    验证并返回有效的 host 地址
    :param host: 用户输入的 host 地址
    :return: 验证后的 host 地址
    """
    api_host: Optional[str] = host
    web_host: Optional[str] = None
    if not host:
        api_host = settings.api_host
        web_host = settings.web_host

    # 此时api_host必然存在
    api_host, web_host, _ = resolve_hosts(api_host=api_host, web_host=web_host)
    assert api_host is not None
    assert web_host is not None

    return api_host, web_host


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

    :c
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

    # 先使用 console 打印提示，因为遮罩读取原生不支持 Rich 的颜色标签渲染
    console.print(prompt_text, end="")

    # 强制刷新输出缓冲区，确保提示语立刻显示
    sys.stdout.flush()

    # 3. 安全读取用户输入（输入时以星号遮罩，避免用户误以为终端卡死）
    with safe.block(message="Failed to read API Key from terminal"):
        try:
            key = utils.prompt_masked()
            return key.strip()
        except (KeyboardInterrupt, EOFError):
            # 优雅处理用户按下 Ctrl+C 或 Ctrl+D 退出的情况，替代旧版的 sys.excepthook
            console.print("\n")  # 换行，防止终端提示符错位
            sys.exit(0)


def welcome(base_url: str, profile: UserProfile):
    """
    登录成功后打印欢迎信息
    :param base_url: 登录地址
    :param profile: 用户档案对象，包含用户信息等数据
    :return:
    """
    username = profile.get("username", "")
    nickname = profile.get("name", "")
    name = nickname or username
    if name:
        console.info(
            "Currently logged in as:",
            Text(name, "yellow"),
            "to",
            Text.assemble((base_url, "green"), ". Use"),
            Text("`swanlab login --relogin`", "bold"),
            "to force relogin",
        )
