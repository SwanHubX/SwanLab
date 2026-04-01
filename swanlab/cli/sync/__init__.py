"""
@author: caddiesnew
@file: __init__.py
@time: 2026/4/1 16:10
@description: CLI Sync 模块：同步本地数据到云端
"""

import click


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
    help="The workspace to sync the logs to. If not specified, uses the default workspace.",
)
@click.option(
    "--project",
    "-p",
    default=None,
    type=str,
    help="The project to sync the logs to. If not specified, uses the default project.",
)
@click.option(
    "--id",
    "-i",
    default=None,
    type=str,
    help="The experiment ID to sync the logs to. Only valid when path is a single directory.",
)
def sync(path, api_key, workspace, project, host, id):
    """Synchronize local logs to the cloud."""
    # TODO: 接入重构后的 sync 逻辑
    pass
