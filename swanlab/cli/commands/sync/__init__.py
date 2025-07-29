"""
@author: cunyue
@file: __init__.py
@time: 2025/6/5 14:03
@description: 同步本地数据到云端
"""

import click

from swanlab.error import KeyFileError
from swanlab.package import HostFormatter, get_key
from swanlab.sync import sync as sync_logs


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
    help="The API key to use for authentication. If not specified, it will use the default API key from the environment."
    "If specified, it will log in using this API key but will not save the key.",
)
@click.option(
    "--host",
    "-h",
    default=None,
    type=str,
    help="The host to sync the logs to. If not specified, it will use the default host.",
)
@click.option(
    "--workspace",
    "-w",
    default=None,
    type=str,
    help="The workspace to sync the logs to. If not specified, it will use the default workspace.",
)
@click.option(
    "--project",
    "-p",
    default=None,
    type=str,
    help="The project to sync the logs to. If not specified, it will use the default project.",
)
@click.option(
    "--id",
    "-i",
    default=None,
    type=str,
    help="The experiment ID to sync the logs to. It can only be used when the path is a single directory."
    "For more details, see https://docs.swanlab.cn/api/cli-swanlab-sync.html",
)
def sync(path, api_key, workspace, project, host, id):
    """
    Synchronize local logs to the cloud.
    """
    # 1. 检查 host 参数是否符合规范，并注入环境变量中
    HostFormatter(host)()
    # 2. 尝试获取 API key
    try:
        if api_key is None:
            api_key = get_key()
    except KeyFileError:
        pass
    # 2. 同步数据
    for p in path:
        sync_logs(p, workspace=workspace, project=project, id=id, api_key=api_key, raise_error=len(path) == 1)
