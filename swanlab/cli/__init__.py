"""
@author: cunyue
@file: __init__.py
@time: 2026/3/5 14:32
@description: SwanLab 命令行交互，在终端运行 `swanlab` 命令进行一些CLI操作
"""

import click

from swanlab.sdk.internal.pkg.helper import get_swanlab_version

from .auth import login, logout, verify
from .converter import convert
from .dashboard import watch
from .mode import disabled, local, offline, online
from .sync import sync


@click.group(invoke_without_command=True)
@click.version_option(get_swanlab_version(), "--version", "-v", message="SwanLab %(version)s")
def cli():
    """SwanLab CLI — Track, compare, and visualize AI experiments."""
    pass


# ---------------------------------- 注册子命令 ----------------------------------
# Auth
# noinspection PyTypeChecker
cli.add_command(login)
# noinspection PyTypeChecker
cli.add_command(logout)
# noinspection PyTypeChecker
cli.add_command(verify)

# Dashboard
# noinspection PyTypeChecker
cli.add_command(watch)

# Converter
# noinspection PyTypeChecker
cli.add_command(convert)

# Sync
# noinspection PyTypeChecker
cli.add_command(sync)

# Mode
# noinspection PyTypeChecker
cli.add_command(offline)
# noinspection PyTypeChecker
cli.add_command(local)
# noinspection PyTypeChecker
cli.add_command(online)
# noinspection PyTypeChecker
cli.add_command(disabled)


__all__ = ["cli"]
