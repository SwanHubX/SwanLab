"""
@author: cunyue
@file: __init__.py
@time: 2025/6/5 14:03
@description: 同步本地数据到云端
"""

import click

from swanlab.api import terminal_login, create_http
from swanlab.sync import sync as sync_logs


@click.command()
@click.argument(
    "path",
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    nargs=1,
    required=True,
)
def sync(path):
    log_info = terminal_login()
    create_http(log_info)
    sync_logs(path, login_required=False)
