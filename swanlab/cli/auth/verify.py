"""
@author: caddiesnew
@file: verify.py
@time: 2026/4/9 20:38
@description: CLI 验证登录状态命令
"""

import sys

import click

from swanlab.sdk.cmd.verify import raw_verify


@click.command()
def verify():
    """Verify the current login status."""
    success = raw_verify()
    if not success:
        sys.exit(1)
