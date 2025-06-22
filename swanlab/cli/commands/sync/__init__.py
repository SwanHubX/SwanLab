"""
@author: cunyue
@file: __init__.py
@time: 2025/6/5 14:03
@description: 同步本地数据到云端
"""

import click

from swanlab.core_python import create_client, auth
from swanlab.error import KeyFileError
from swanlab.package import get_key, HostFormatter
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
def sync(path, api_key, workspace, project, host):
    """
    Synchronize local logs to the cloud.
    """
    # 1. 创建 http 对象
    # 1.1 检查host是否合法，并格式化，注入到环境变量中
    HostFormatter(host)()
    # 1.2 如果输入了 api-key， 使用此 api-key 登录但不保存数据
    try:
        api_key = get_key() if api_key is None else api_key
    except KeyFileError:
        pass
    for p in path:
        # 1.3 登录，创建 http 对象
        log_info = auth.terminal_login(api_key=api_key, save_key=False)
        create_client(log_info)
        # 2. 同步日志
        sync_logs(p, workspace=workspace, project_name=project, login_required=False, raise_error=len(path) == 1)
