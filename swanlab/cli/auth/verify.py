"""
@author: caddiesnew
@file: verify.py
@time: 2026/4/9 20:38
@description: CLI 验证登录状态命令
"""

import sys

import click
from rich.text import Text

from swanlab import sdk


@click.command()
@click.option("--local", is_flag=True, help="Verify local login status (check .swanlab in current directory)")
def verify(local: bool):
    """Verify the current login status."""

    # 1. 检查本地是否存在 API Key
    nrc_path = sdk.cmd_utils.get_nrc_path(save="local" if local else "root")
    if not nrc_path.exists():
        sdk.pkg.console.error("You are not logged in. Use `swanlab login` to login.")
        sys.exit(1)

    # 2. 读取 API Key 并验证
    netrc_result = sdk.pkg.nrc.read(nrc_path)
    if not netrc_result:
        sdk.pkg.console.error("You are not logged in. Use `swanlab login` to login.")
        sys.exit(1)
    api_key, api_host, web_host = netrc_result
    login_resp = sdk.pkg.client.login_by_api_key(base_url=api_host + "/api", api_key=api_key)

    if login_resp is None:
        sdk.pkg.console.error(
            "Verification failed. Please check if your API key is correct and try again. Use `swanlab login` to login if you haven't done so.",
        )
        sys.exit(1)

    # 3. 展示验证结果
    username = login_resp.get("userInfo", {}).get("username", "unknown")
    sdk.pkg.console.info("You are logged into", Text(web_host, style="blue"), "as", Text(username, "bold"), sep=" ")
