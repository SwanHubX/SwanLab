"""
@author: cunyue
@file: login.py
@time: 2026/3/6 22:24
@description: swanlab.login 方法，登录到 SwanLab 平台
"""

from typing import Optional

from swanlab.exceptions import AuthenticationError
from swanlab.sdk.internal import apikey
from swanlab.sdk.internal.context import RunConfig, RunContext, get_context, has_context, use_temp_context
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg.scope import Scope
from swanlab.sdk.internal.settings import Settings, settings
from swanlab.sdk.pkg import console, helper
from swanlab.sdk.typings.core_python.api.bootstrap import LoginResponse


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

    [Note that] this function should be called before `swanlab.init`, like swanlab.merge_settings.

    :param api_key: str, optional
        Your SwanLab authentication key. If not provided, the SDK will attempt to
        read it from the local environment or prompt you to input it interactively.
    :param relogin: bool, optional
        If True, forces a re-authentication and overwrites the existing API key.
        When `relogin` is True and you already logged in, this function will return immediately.
        Defaults to False.
    :param host: str, optional
        The API host URL. If not provided, the default SwanLab cloud host will be used.
        This parameter equivalent to setting the `SWANLAB_API_HOST` environment variable.
        If you set this parameter, settings and environment variables will be ignored and updated.
    :param save: bool, optional
        Whether to save the provided api_key to the local netrc/credential file
        for future sessions. Defaults to False.
    :param timeout: int, optional
        Timeout in seconds for the login network request. Defaults to 10.

    :raises RuntimeError: If the environment does not support interactive input.
    :raises AuthenticationError: If login fails.
    :return: Returns True if login was successful, False otherwise.
    """
    # 1. 如果已经登录且不需要重新登录，则直接返回
    if not relogin and client.exists():
        console.info("Already logged in, Skipping login. If you want to relogin, use `swanlab.login(relogin=True)`")
        return True
    # 如果已经初始化了运行上下文，则不允许重新登录
    if has_context():
        console.error("Cannot relogin while SwanLab Context is active. Please clear the context first.")
        return False
    # 2. 获取当前配置
    api_key = api_key or settings.api_key
    api_host = host or settings.api_host
    web_host = host or settings.web_host
    login_settings = Settings.model_validate({"api_key": api_key, "api_host": api_host, "web_host": web_host})
    # 临时使用运行上下文，在结束后清除
    with use_temp_context(RunContext(config=RunConfig(settings=login_settings))) as ctx:
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


@helper.rich.with_loading_animation()
def create_client(timeout: int = 10):
    ctx = get_context()
    assert ctx.config.settings.api_key is not None, "API Key not provided"
    return client.new(ctx.config.settings.api_key, ctx.config.settings.api_host, timeout=timeout)
