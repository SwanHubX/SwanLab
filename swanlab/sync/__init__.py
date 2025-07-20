import os.path
from sys import stdout
from typing import Literal, Union

from rich.status import Status

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from ..core_python import create_client, get_client
from ..core_python.auth.providers.api_key import code_login
from ..data.porter import DataPorter, Mounter
from ..formatter import check_proj_name_format, check_run_id_format
from ..log import swanlog


def sync(
    dir_path: str,
    workspace: str = None,
    project: str = None,
    id: Union[str, Literal['auto', 'new']] = None,
    api_key: str = None,
    raise_error: bool = True,
):
    """
    Syncs backup files to the cloud. Before syncing, you must log in.
    :param dir_path: The directory path to sync.
    :param id: The ID of the backup to sync. Use cases:
        - None: Equal to 'new'
        - new: Create a new experiment with a new ID.
        - auto: Create (Resume) the experiment with the ID from the backup file.
        - str: Use the specified ID to sync the logs.
    :param workspace: The workspace to sync the logs to. If not specified, it will use the default workspace.
    :param project: The project to sync the logs to. If not specified, it will use the default project.
    :param raise_error: Whether to raise an error if error occurs when syncing.
    :param api_key: If provided, swanlab will sync using this API key. Or you need login first before run this function.
    """
    # 0. 参数检查
    # 0.1 检查项目名称
    project and check_proj_name_format(project)
    # 0.2 检查实验 ID
    if id is None:
        id = "new"
    if id not in ['new', 'auto']:
        check_run_id_format(id)
    # 1. 根据 api key 登录
    try:
        client = get_client()
    except ValueError:
        client = None
    # api key 存在，则尝试创建客户端
    if api_key:
        client = create_client(code_login(api_key, save_key=False))
    # 2. 开始同步
    try:
        assert os.path.exists(dir_path), f"Directory {dir_path} does not exist."
        stdout.flush()
        with Status("🔁 Syncing...", spinner="dots"):
            with DataPorter().open_for_sync(run_dir=dir_path) as porter:
                proj, exp = porter.parse()
                assert client is not None, "Please log in first before using sync."
                with Mounter() as mounter:
                    run_store = mounter.run_store
                    # 设置项目信息
                    run_store.project = project or proj.name
                    run_store.workspace = workspace or proj.workspace
                    run_store.visibility = proj.public
                    # 设置实验信息
                    run_store.run_name = exp.name
                    run_store.description = exp.description
                    run_store.tags = exp.tags
                    run_store.run_colors = exp.colors
                    # 设置实验 id 和 resume 模式
                    # a. id 为 new，则 resume 为 never, id 为 None
                    # b. id 为 auto，则 resume 为 allow, id 为 exp.id
                    # c. 其他情况，则 resume 为 must, id 为 id
                    if id == "new":
                        run_store.resume = 'never'
                        run_store.run_id = None
                    elif id == "auto":
                        run_store.resume = 'allow'
                        run_store.run_id = exp.id
                    else:
                        run_store.resume = 'must'
                        run_store.run_id = id
                    # 创建实验会话
                    mounter.execute()
                    # 同步
                    _ = porter.synchronize()
        swanlog.info("🚀 Sync completed, View run at ", client.web_exp_url)
    except Exception as e:
        if raise_error:
            raise e
        else:
            return swanlog.error(f"😢 Error syncing: {e}")
