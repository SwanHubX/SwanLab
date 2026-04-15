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
    success = sdk.logout_cli(force=force, save="local" if local else "root")
    if not success:
        sys.exit(1)
