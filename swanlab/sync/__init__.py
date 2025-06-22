import os.path
from sys import stdout

from rich.status import Status

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from ..core_python import get_client
from ..data.namer import generate_colors
from ..data.porter import DataPorter
from ..log import swanlog


def sync(
    dir_path: str,
    workspace: str = None,
    project_name: str = None,
    raise_error: bool = True,
    login_required: bool = True,
):
    """
    Syncs backup files to the cloud. Before syncing, you must log in.
    :param dir_path: The directory path to sync.
    :param workspace: The workspace to sync the logs to. If not specified, it will use the default workspace.
    :param project_name: The project to sync the logs to. If not specified, it will use the default project.
    :param raise_error: Whether to raise an error if error occurs when syncing.
    :param login_required: Whether login is required before syncing, just for debugging.
    """
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
                project, experiment = porter.parse()
                assert client is not None, "Please log in first, use `swanlab login` to log in."
                client.mount_project(
                    name=project_name or project.name,
                    username=workspace or project.workspace,
                    public=project.public,
                )
                colors = generate_colors(client.history_exp_count)
                client.mount_exp(
                    exp_name=experiment.name,
                    colors=colors,
                    description=experiment.description,
                    tags=experiment.tags,
                )
                success = porter.synchronize()
                # 3.5 更新实验状态
                client.update_state(success=success)
        swanlog.info("🚀 Sync completed, View run at ", client.web_exp_url)
    except Exception as e:
        if raise_error:
            raise e
        else:
            return swanlog.error(f"😢 Error syncing: {e}")
