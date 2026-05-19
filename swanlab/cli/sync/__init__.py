"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/1 16:10
@description: CLI Sync 模块：同步本地数据到云端
"""

from pathlib import Path
from typing import Optional, Tuple

import click

from swanlab.sdk import Settings, pkg
from swanlab.sdk import sync as sync_cmd


@click.command()
@click.argument(
    "path",
    nargs=-1,
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    required=True,
)
@click.option(
    "--api-key",
    "-k",
    default=None,
    type=str,
    help="The API key to use for authentication. If not specified, uses the default API key.",
)
@click.option(
    "--host",
    "-h",
    default=None,
    type=str,
    help="The host to sync the logs to. If not specified, uses the default host.",
)
@click.option(
    "--workspace",
    "-w",
    default=None,
    type=str,
    help="The workspace to sync the logs to. If not specified, uses the previous workspace.",
)
@click.option(
    "--project",
    "-p",
    default=None,
    type=str,
    help="The project to sync the logs to. If not specified, uses the previous project.",
)
@click.option(
    "--id",
    "-i",
    default=None,
    type=str,
    help="The run id to sync the logs to.",
)
def sync(
    path: Tuple[str],
    api_key: Optional[str],
    workspace: Optional[str],
    project: Optional[str],
    host: Optional[str],
    id: Optional[str],
):
    """Synchronize local logs to the cloud."""
    model = pkg.helper.strip_none(
        {
            "api_key": api_key,
            "api_host": host,
            "project": {"workspace": workspace, "name": project},
            "run": {"id": id},
        },
        strip_empty_str=True,
        strip_empty_dict=True,
    )
    settings = Settings.model_validate(model)
    for p in path:
        sync_cmd(Path(p), settings=settings)
