"""
@author: cunyue
@file: login.py
@time: 2026/3/6 22:24
@description: swanlab.login 方法，登录到 SwanLab 平台
"""

from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.cmd.helper import with_cmd_lock
from swanlab.sdk.internal import apikey
from swanlab.sdk.internal.context import RunConfig, RunContext, get_context, has_context, use_temp_context
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.scope import Scope
from swanlab.sdk.internal.run import has_run
from swanlab.sdk.internal.settings import Settings, settings
from swanlab.sdk.typings.core_python.api.bootstrap import LoginResponse
from swanlab.sdk.utils import helper


@with_cmd_lock
def login(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: bool = False,
    timeout: int = 10,
) -> bool:
    """
    Login to SwanLab Cloud.

    This function authenticates your environment with SwanLab. If an API key is
    already configured locally and `relogin` is False, this function will do nothing.

    [Note that] this function should be called before `swanlab.init`, like `swanlab.merge_settings`.

    Examples:
    ---------
    >>> import swanlab
    >>> # 1. Basic login with explicit API key
    >>> swanlab.login("your_api_key")
    >>>
    >>> # 2. Interactive login (prompts for input if no key is found in env/config)
    >>> swanlab.login()
    >>>
    >>> # 3. Force re-login and save the new credential to local .netrc
    >>> swanlab.login("new_api_key", relogin=True, save=True)
    >>>
    >>> # 4. Login to a private SwanLab deployment
    >>> swanlab.login("private_key", host="https://private-swanlab.com")
    >>>
    >>> # Now you can start your experiment
    >>> swanlab.init()

    :param api_key: Your SwanLab authentication key. If not provided, the SDK will attempt to
        read it from the local environment or prompt you to input it interactively.

    :param relogin: If True, forces a re-authentication and overwrites the existing API key.
        When `relogin` is False and you are already logged in, this function will return immediately.
        Defaults to False.

    :param host: The API host URL. If not provided, the default SwanLab cloud host will be used.
        This parameter is equivalent to setting the `SWANLAB_API_HOST` environment variable.
        If you set this parameter, global settings and environment variables will be overwritten.

    :param save: Whether to save the provided api_key to the local netrc/credential file
        for future sessions. Defaults to False.

    :param timeout: Timeout in seconds for the login network request. Defaults to 10.

    :raises RuntimeError: If the environment does not support interactive input and no key is provided.

    :raises AuthenticationError: If the login request fails (e.g., invalid key or network issue).

    :return: Returns True if login was successful, False otherwise.
    """
    if has_run():
        console.error("Cannot login while SwanLab Run is active. Please finish the run first.")
        return False
    # 1. 判断是是否允许重新登录
    # 如果已经初始化了运行上下文，则不允许重新登录
    if has_context():
        console.error("Cannot relogin while SwanLab Context is active. Please clear the context first.")
        return False
    # 如果已经登录且不需要重新登录，则直接返回
    if not relogin and client.exists():
        console.info("Already logged in, Skipping login. If you want to relogin, use `swanlab.login(relogin=True)`")
        return True
    elif client.exists() and relogin:
        client.reset()
    # 2. 获取当前配置
    # 如果提供了host且没有http前缀添加https://前缀
    if host is not None:
        host = host.strip().rstrip("/")
        if not host.startswith(("http://", "https://")):
            host = f"https://{host}"
        parsed = urlparse(host)
        host = f"{parsed.scheme}://{parsed.netloc}"
    # 与settings中同步
    api_key = api_key or settings.api_key
    # 如果输入了host，检查与settings中的api_host是否一致，不一致则清空api_key
    # 额外还需要判断是否存在本地.netrc文件，如果存在则不覆盖
    if host and api_key and settings.api_host != host and apikey.exists_locally():
        api_key = None
    api_host = host or settings.api_host
    # 防止用户通过 login(host="api.swanlab.cn") 强行覆盖 web_host
    if host and "api.swanlab.cn" in host:
        web_host = settings.web_host  # 保持官方默认 web_host
    else:
        web_host = host or settings.web_host
    login_settings = Settings.model_validate({"api_key": api_key, "api_host": api_host, "web_host": web_host})
    # 临时使用运行上下文，在结束后清除
    fake_run_dir = Path.cwd()
    with use_temp_context(RunContext(config=RunConfig(settings=login_settings, run_dir=fake_run_dir))) as ctx:
        # 如果 API Key 不存在，则提示用户输入
        if api_key is None:
            api_key = apikey.prompt()
        # 3. 进入登录流程
        ctx.config.settings.merge_settings({"api_key": api_key})
        with Scope() as scope:
            create_client(timeout=timeout)
            assert client.exists(), "Failed to create client"
            login_resp: Optional[LoginResponse] = scope.get("login_resp", None)
            if login_resp is None:
                raise AuthenticationError("Failed to login, please check your API Key or network connection.")
            if save:
                apikey.save(username=web_host, api_key=api_key, host=ctx.config.settings.api_host)
        # 4. 将登录设置合并到全局配置中
        settings.merge_settings(ctx.config.settings)
        return True


def interactive_login(
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
    try:
        # 首次尝试登录，复用原子接口
        return login(api_key=api_key, relogin=relogin, host=host, save=save, timeout=timeout)
    except AuthenticationError as e:
        # 如果全局配置禁用了交互模式，直接抛出异常
        if not settings.interactive:
            raise e
        console.error(str(e))

    # 进入容错重试循环
    while True:
        try:
            # 重新要求用户输入新的 Key
            new_key = apikey.prompt()
            # 必须设置 relogin=True，确保重置底层的旧 client 状态
            return login(api_key=new_key, relogin=True, host=host, save=save, timeout=timeout)
        except AuthenticationError as e:
            console.error(str(e))
        except (KeyboardInterrupt, EOFError):
            console.info("\nLogin cancelled by user.")
            return False


@helper.rich.with_loading_animation()
def create_client(timeout: int = 10):
    ctx = get_context()
    assert ctx.config.settings.api_key is not None, "API Key not provided"
    return client.new(ctx.config.settings.api_key, ctx.config.settings.api_host, timeout=timeout)
