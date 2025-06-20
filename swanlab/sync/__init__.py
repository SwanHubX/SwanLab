import os.path
from sys import stdout

from rich.status import Status

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from ..core_python import get_client, uploader
from ..data.namer import generate_colors
from ..log import swanlog
from swanlab.data.backup import BackupHandler, DataStore
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
            # 0. 检查备份文件
            assert (
                header.backup_type == "DEFAULT"
            ), f"Backup type mismatch: {header.backup_type}, please update your swanlab package."
            # 1. 聚合信息
            # 1.1 聚合日志
            log_model_dict = {"INFO": [], "WARN": [], "ERROR": []}
            for log in logs:
                log_model = log.to_log_model()
                log_model_dict[log_model['level']].append(log_model)
            # 1.2 生成运行时
            runtime_model = runtime.to_file_model(os.path.join(dir_path, "files"))
            # 1.3 集合列
            column_models = [c.to_column_model() for c in columns]
            # 1.4 集合指标
            scalar_models = [scalar.to_scalar_model() for scalar in scalars]
            # 1.5 集合媒体
            media_models = [media.to_media_model(os.path.join(dir_path, "media")) for media in medias]
            assert client is not None, "Please log in first, use `swanlab login` to log in."
    except Exception as e:
        if raise_error:
            raise e
        else:
            return swanlog.error(f"❌  Error parsing backup file: {e}")
    # 第二部分，上传数据到云端
    try:
        # 3. 上传数据
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
            # 3.2 上传日志
            for key in log_model_dict:
                if len(log_model_dict[key]) > 0:
                    uploader.upload_logs(log_model_dict[key])
            # 3.3 上传运行时
            uploader.upload_files([runtime_model])
            # 3.4 上传列、标量、媒体
            uploader.upload_columns(column_models)
            uploader.upload_scalar_metrics(scalar_models)
            uploader.upload_media_metrics(media_models)
            # 3.5 更新实验状态
            client.update_state(success=footer.success if footer else False)
    except Exception as e:
        if raise_error:
            raise e
        else:
            swanlog.error(f"❌  Error uploading data: {e}")
            return
    swanlog.info("🚀 Sync completed, View run at ", client.web_exp_url)
