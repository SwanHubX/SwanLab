"""
@author: caddiesnew
@file: logout.py
@time: 2026/4/9 20:38
@description: CLI 登出命令
"""

import sys

import click

from swanlab import sdk


@click.command()
@click.option(
    "--force",
    "-f",
    is_flag=True,
    default=False,
    help="Force logout without confirmation prompt.",
)
@click.option("--local", is_flag=True, help="Logout from local login (remove .swanlab in current directory)")
def logout(force: bool, local: bool):
    """Logout from the SwanLab cloud."""

    # 1. 检查本地是否存在 API Key
    nrc_path = sdk.cmd_utils.get_nrc_path(save="local" if local else "root")
    if not nrc_path.exists():
        sdk.pkg.console.error("You are not logged in. Use `swanlab login` to login.")
        sys.exit(1)

    # 2. 交互式确认（非 force 模式下）
    if not force:
        if sdk.settings.interactive:
            sdk.pkg.console.info("Are you sure you want to logout? (y/N): ", end="")
            try:
                confirm = input()
            except (KeyboardInterrupt, EOFError):
                sdk.pkg.console.info("\nLogout canceled.")
                sys.exit(1)
            if confirm.lower() != "y":
                sdk.pkg.console.info("Logout canceled.")
                sys.exit(1)
        else:
            sdk.pkg.console.info("Use 'swanlab logout --force' to skip confirmation in non-interactive mode.")
            sys.exit(1)
    sdk.pkg.nrc.remove(nrc_path)
    sdk.pkg.console.info("Logout successfully. You can use `swanlab login` to login again.")
