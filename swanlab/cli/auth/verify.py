"""
@author: caddiesnew
@file: verify.py
@time: 2026/4/9 20:38
@description: CLI 验证登录状态命令
"""

import sys

import click

from swanlab.sdk.cmd.verify import verify_raw


@click.command()
def verify():
    """Verify the current login status."""
    success = verify_raw()
    if not success:
        sys.exit(1)
