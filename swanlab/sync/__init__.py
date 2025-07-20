import os.path
from sys import stdout
from typing import Literal, Union

from rich.status import Status

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from ..core_python import get_client
from ..data.porter import DataPorter, Mounter
from ..formatter import check_proj_name_format
from ..log import swanlog


def sync(
    dir_path: str,
    workspace: str = None,
    project: str = None,
    id: Union[str, Literal['auto']] = None,
    raise_error: bool = True,
    login_required: bool = True,
):
    """
    Syncs backup files to the cloud. Before syncing, you must log in.
    :param dir_path: The directory path to sync.
    :param id: The ID of the backup to sync. Use cases:
        - None(default): Create a new experiment with a new ID.
        - 'auto': Create a new experiment with the ID from the backup file.
        - str: Use the specified ID to sync the logs.
    :param workspace: The workspace to sync the logs to. If not specified, it will use the default workspace.
    :param project: The project to sync the logs to. If not specified, it will use the default project.
    :param raise_error: Whether to raise an error if error occurs when syncing.
    :param login_required: Whether login is required before syncing, just for debugging.
    """
    # 检查格式
    project and check_proj_name_format(project)
    if id != 'auto':
        id and check_proj_name_format(id)
    # 第一部分，处理备份文件，读取备份信息到内存中
    try:
        assert os.path.exists(dir_path), f"Directory {dir_path} does not exist."
        try:
            client = get_client()
        except ValueError:
            client = None
            assert not login_required, "Please log in first, use `swanlab login` to log in."
        stdout.flush()
        with Status("🔁 Syncing...", spinner="dots"):
            with DataPorter().open_for_sync(run_dir=dir_path) as porter:
                proj, exp = porter.parse()
                assert client is not None, "Please log in first, use `swanlab login` to log in."
                with Mounter() as mounter:
                    run_store = mounter.run_store
                    # 设置项目信息
                    run_store.project = project or proj.name
                    run_store.workspace = workspace or proj.workspace
                    run_store.visibility = proj.public
                    # 设置实验信息
                    run_store.run_id = exp.id
                    run_store.run_name = exp.name
                    run_store.description = exp.description
                    run_store.tags = exp.tags
                    run_store.run_colors = exp.colors
                    # 创建实验会话
                    mounter.execute()
                    # 同步
                    success = porter.synchronize()
                    # 更新实验状态
                    client.update_state(success=success)
        swanlog.info("🚀 Sync completed, View run at ", client.web_exp_url)
    except Exception as e:
        if raise_error:
            raise e
        else:
            return swanlog.error(f"😢 Error syncing: {e}")
