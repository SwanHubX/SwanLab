"""
@author: cunyue
@file: __init__.py
@time: 2025/6/5 14:03
@description: 同步本地数据到云端
"""

import click

from swanlab.core_python import create_client, auth
from swanlab.error import KeyFileError
from swanlab.formatter import check_proj_name_format, check_run_id_format
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
@click.option(
    "--id",
    "-i",
    default=None,
    type=str,
    help="The experiment ID to sync the logs to. If not specified, it will use the default experiment ID. "
    "It can only be used when the path is a single directory. ",
)
@click.option(
    "--resume",
    "-r",
    is_flag=True,
    default=False,
    help="If set, it will resume the run if it exists, otherwise it will create a new run. "
    "Resume can not be used with --id, --workspace, --project options, and it can only be used when the path is a single directory.",
)
def sync(path, api_key, workspace, project, host, id, resume):
    """
    Synchronize local logs to the cloud.
    """
    # 0. 参数检查
    # 0.1 检查 project, id 参数是否符合规范
    project and check_proj_name_format(project)
    id and check_run_id_format(id)
    # 0.2 如果 resume 被设置，则不允许使用 id, workspace, project 参数
    if resume and (id or workspace or project):
        raise click.UsageError("The --resume option cannot be used with --id, --workspace, or --project options.")
    # 0.3 如果 id 或者 resume 被设置，则 path 必须是单个目录
    if (id or resume) and len(path) != 1:
        raise click.UsageError("The --id or --resume option can only be used with a single directory path.")
    id = "auto" if resume else id
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
        # 在后端，每一个 http 对象对应一个实验会话，所以如果有多个路径，则需要为每个路径创建一个 http 对象
        log_info = auth.terminal_login(api_key=api_key, save_key=False)
        create_client(log_info)
        # 2. 同步日志
        sync_logs(p, workspace=workspace, project=project, id=id, login_required=False, raise_error=len(path) == 1)
