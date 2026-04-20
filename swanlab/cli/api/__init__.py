"""
@author: Nexisato
@file: __init__.py
@time: 2026/4/20
@description: CLI API 子命令 — 通过命令行调用 SwanLab 公共查询 API
"""

import click

from .experiment import get_run
from .project import get_project
from .user import get_user
from .workspace import get_workspace


@click.group("api")
def api_cli():
    """Generic SwanLab API requests."""
    pass


api_cli.add_command(get_project)
api_cli.add_command(get_run)
api_cli.add_command(get_workspace)
api_cli.add_command(get_user)


__all__ = ["api_cli"]
