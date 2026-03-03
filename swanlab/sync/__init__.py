import os.path
from sys import stdout
from typing import Literal, Union

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from .sync_utils import set_run_store, SyncProgress
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
        Attention: If you not provide api-key, you need login every time you run this function.
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
        with SyncProgress() as pbar:
            with DataPorter().open_for_sync(run_dir=dir_path, progress_callback=pbar.update) as porter:
                proj, exp = porter.parse()
                assert client is not None, "Please log in first before using sync."
                with Mounter() as mounter:
                    run_store = mounter.run_store
                    # 设置运行存储相关信息
                    set_run_store(run_store, proj, exp, project, workspace, id)
                    # 创建实验会话
                    mounter.execute()
                    # 显示最终 URL（在上传进度条开始前）
                    swanlog.info(f"🔁 Syncing data to cloud... (View at {client.web_exp_url})")
                    # 同步
                    porter.synchronize(set_total_callback=pbar.set_total)
        swanlog.info("🚀 Sync completed, View run at ", client.web_exp_url)
    except Exception as e:
        if raise_error:
            raise e
        else:
            return swanlog.error(f"😢 Error syncing: {e}")
