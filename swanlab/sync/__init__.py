import os.path
from sys import stdout

from rich.status import Status

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from ..core_python import get_client
from ..data.namer import generate_colors
from swanlab.transfers import ProtoV0Transfer
from ..log import swanlog
from swanlab.log.backup import BackupHandler, DataStore
from swanlab.proto.v0 import ModelsParser


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
        with Status("🛠️Parsing...", spinner="dots"):
            assert os.path.exists(dir_path), f"Directory {dir_path} does not exist."
            file_path = os.path.join(dir_path, BackupHandler.BACKUP_FILE)
            assert os.path.exists(file_path), f"Can not find backup file {BackupHandler.BACKUP_FILE} in {dir_path}."
            try:
                client = get_client()
            except ValueError:
                client = None
                assert not login_required, "Please log in first, use `swanlab login` to log in."
            stdout.flush()
            ds = DataStore()
            ds.open_for_scan(file_path)
            # 根据类型分类日志文件
            with ModelsParser() as models_parser:
                for record in ds:
                    if record is None:
                        continue
                    models_parser.parse_record(record)
            header, project, experiment, logs, runtime, columns, scalars, medias, footer = models_parser.get_parsed()
            # 检查备份文件
            assert (
                header.backup_type == "DEFAULT"
            ), f"Backup type mismatch: {header.backup_type}, please update your swanlab package."
            assert client is not None, "Please log in first, use `swanlab login` to log in."
    except Exception as e:
        if raise_error:
            raise e
        else:
            return swanlog.error(f"❌  Error parsing backup file: {e}")
    # 第二部分，上传数据到云端
    try:
        # 3. 上传数据
        transfer = ProtoV0Transfer(media_dir=os.path.join(dir_path, "media"), file_dir=os.path.join(dir_path, "files"))
        # 3.1 创建项目与实验
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
        with Status("🔁 Syncing...", spinner="dots"):
            transfer.publish_file(runtime)
            [transfer.publish_log(log) for log in logs]
            [transfer.publish_column(column) for column in columns]
            [transfer.publish_scalar(scalar) for scalar in scalars]
            [transfer.publish_media(media) for media in medias]
            transfer.join()
            # 3.5 更新实验状态
            client.update_state(success=footer.success if footer else False)
    except Exception as e:
        if raise_error:
            raise e
        else:
            swanlog.error(f"❌  Error uploading data: {e}")
            return
    swanlog.info("🚀 Sync completed, View run at ", client.web_exp_url)
