"""
@author: caddiesnew
@file: verify.py
@time: 2026/4/9 21:00
@description: verify 方法，验证当前登录状态
"""

from rich.text import Text

from swanlab.sdk.cmd import utils
from swanlab.sdk.internal.pkg import console, nrc
from swanlab.sdk.internal.pkg.client.bootstrap import login_by_api_key
from swanlab.sdk.internal.settings import settings
from swanlab.sdk.typings.cmd import LoginType

__all__ = ["verify_cli"]


def verify_cli(save: LoginType = "root") -> bool:
    # 1. 检查本地是否存在 API Key
    nrc_path = utils.get_nrc_path(save)
    if not nrc_path.exists():
        console.info("You are not logged in. Use `swanlab login` to login.")
        return False

    # 2. 读取 API Key 并验证
    netrc_result = nrc.read(nrc_path)
    if not netrc_result:
        console.error("You are not logged in. Use `swanlab login` to login.")
        return False
    api_key, api_host, _ = netrc_result
    login_resp = login_by_api_key(base_url=api_host + "/api", api_key=api_key)

    if login_resp is None:
        console.error(
            "Verification failed. Please check if your API key is correct and try again. Use `swanlab login` to login if you haven't done so.",
        )
        return False

    # 3. 展示验证结果
    username = login_resp.get("userInfo", {}).get("username", "unknown")
    console.info("You are logged into", Text(settings.web_host, style="blue"), "as", Text(username, "bold"), sep=" ")
    return True
