"""
@author: caddiesnew
@file: verify.py
@time: 2026/4/9 21:00
@description: swanlab.verify 方法，验证当前登录状态
"""

from rich.text import Text

from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal import apikey
from swanlab.sdk.internal.core_python.api.bootstrap import login_by_api_key
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.settings import settings


@with_cmd_lock
@without_run("verify")
def verify() -> bool:
    """Verify the current login status.

    This function checks whether the locally stored API key is still valid
    by attempting to authenticate with the SwanLab server.

    :return: True if the login is valid, False otherwise.

    :raises RuntimeError: If called while a run is active.

    Examples:

        Verify current login status:

        >>> import swanlab
        >>> swanlab.verify()
    """
    return raw_verify()


def raw_verify() -> bool:
    # 1. 检查本地是否存在 API Key
    if not apikey.exists():
        console.info("You are not logged in. Use `swanlab login` to login.")
        return False

    # 2. 读取 API Key 并验证
    try:
        key = apikey.get()
    except FileNotFoundError:
        console.info("You are not logged in. Use `swanlab login` to login.")
        return False

    # 3. 调用底层鉴权接口验证 Key 有效性
    api_host = settings.api_host
    base_url = api_host.rstrip("/") + "/api" if not api_host.endswith("/api") else api_host
    login_resp = login_by_api_key(base_url=base_url, api_key=key)

    if login_resp is None:
        console.info("Your API key is invalid or expired. Use `swanlab login --relogin` to login again.")
        return False

    # 4. 展示验证结果
    username = login_resp.get("userInfo", {}).get("username", "unknown")
    console.info("You are logged into", Text(settings.web_host, style="blue"), "as", Text(username, "bold"), sep=" ")
    return True
