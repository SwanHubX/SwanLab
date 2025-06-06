import os.path

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from ..api import get_http
from ..api.upload import upload_logs, upload_files, upload_columns, upload_scalar_metrics, upload_media_metrics
from ..data.namer import generate_colors
from ..log.backup import BackupHandler
from ..log.backup.datastore import DataStore
from ..log.backup.models import ModelsParser, Runtime


def sync(dir_path: str, login_required: bool = True):
    """
    Syncs backup files to the cloud. Before syncing, you must log in.
    :param dir_path: The directory path to sync.
    :param login_required: Whether login is required before syncing, just for debugging.
    """
    assert os.path.exists(dir_path), f"Directory {dir_path} does not exist."
    file_path = os.path.join(dir_path, BackupHandler.BACKUP_FILE)
    assert os.path.exists(file_path), f"Can not find backup file {BackupHandler.BACKUP_FILE} in {dir_path}."
    try:
        http = get_http()
    except ValueError:
        http = None
        assert not login_required, "Please log in first, use `swanlab login` to log in."
    ds = DataStore()
    ds.open_for_scan(file_path)
    # 根据类型分类日志文件
    with ModelsParser() as models_parser:
        for record in ds:
            if record is None:
                continue
            models_parser.parse_record(record)
    header, project, experiment, logs, runtimes, columns, scalars, medias, footer = models_parser.get_parsed()
    # 0. 检查备份文件版本
    assert (
        header.version == BackupHandler.BACKUP_VERSION
    ), f"Backup file version mismatch: {header.version}, please update your swanlab package."
    assert (
        header.backup_type == "DEFAULT"
    ), f"Backup type mismatch: {header.backup_type}, please update your swanlab package."
    # 1. 聚合信息
    # 1.1 聚合日志
    log_model_dict = {"INFO": [], "WARN": [], "ERROR": []}
    for log in logs:
        log_model = log.to_log_model()
        log_model_dict[log_model['level']].append(log_model)
    # 1.2 集合运行时
    runtime = Runtime(
        conda_filename=None,
        requirements_filename=None,
        metadata_filename=None,
        config_filename=None,
    )
    for r in runtimes:
        runtime.conda_filename = r.conda_filename
        runtime.requirements_filename = r.requirements_filename
        runtime.metadata_filename = r.metadata_filename
        runtime.config_filename = r.config_filename
    runtime_model = runtime.to_file_model(os.path.join(dir_path, "files"))
    # 1.3 集合列
    column_models = [c.to_column_model() for c in columns]
    # 1.4 集合指标
    scalar_models = [scalar.to_scalar_model() for scalar in scalars]
    # 1.5 集合媒体
    media_models = [media.to_media_model(os.path.join(dir_path, "media")) for media in medias]

    assert http is not None, "Please log in first, use `swanlab login` to log in."
    # 3. 上传数据
    # 3.1 创建项目与实验
    http.mount_project(name=project.name, username=project.workspace, public=project.public)
    colors = generate_colors(http.history_exp_count)
    http.mount_exp(
        exp_name=experiment.name,
        colors=colors,
        description=experiment.description,
        tags=experiment.tags,
    )
    # 3.2 上传日志
    for key in log_model_dict:
        if len(log_model_dict[key]) > 0:
            upload_logs(log_model_dict[key])
    # 3.3 上传运行时
    upload_files([runtime_model])
    # 3.4 上传列、标量、媒体
    upload_columns(column_models)
    upload_scalar_metrics(scalar_models)
    upload_media_metrics(media_models)
    # 3.5 更新实验状态
    http.update_state(success=footer.success if footer else False)
