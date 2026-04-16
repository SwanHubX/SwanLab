"""
@author: caddiesnew
@file: verify.py
@time: 2026/4/9 20:38
@description: CLI 验证登录状态命令
"""

import sys

import click

from swanlab import sdk


@click.command()
@click.option("--local", is_flag=True, help="Verify local login status (check .swanlab in current directory)")
def verify(local: bool):
    """Verify the current login status."""
    success = sdk.verify_cli(save="local" if local else "root")
    if not success:
        sys.exit(1)
