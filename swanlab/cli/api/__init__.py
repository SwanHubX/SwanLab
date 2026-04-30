"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/20
@description: CLI API 子命令 — 通过命令行调用 SwanLab 公共查询 API
"""

import click

from .experiment import run_cli
from .project import project_cli
from .selfhosted import selfhosted_cli
from .user import user_cli
from .workspace import workspace_cli


@click.group("api")
def api_cli():
    """Generic SwanLab API requests with CLI."""
    pass


api_cli.add_command(project_cli)
api_cli.add_command(run_cli)
api_cli.add_command(workspace_cli)
api_cli.add_command(selfhosted_cli)
api_cli.add_command(user_cli)


__all__ = ["api_cli"]
