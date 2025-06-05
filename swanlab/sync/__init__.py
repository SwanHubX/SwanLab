import os.path

from .mlflow import sync_mlflow
from .tensorboard import sync_tensorboardX, sync_tensorboard_torch
from .wandb import sync_wandb

__all__ = ["sync_wandb", "sync_tensorboardX", "sync_tensorboard_torch", "sync_mlflow", "sync"]

from ..api import get_http
from ..log.backup import BackupHandler
from ..log.backup.datastore import DataStore


def sync(dir_path: str, login_required: bool = True):
    """
    Syncs backup files to the cloud. Before syncing, you must log in.
    :param dir_path: The directory path to sync.
    :param login_required: Whether login is required before syncing, just for debugging.
    """
    assert os.path.exists(dir_path), f"Directory {dir_path} does not exist."
    file_path = os.path.join(dir_path, BackupHandler.BACKUP_FILE)
    assert os.path.exists(file_path), f"Not a valid swanlab sync directory: {dir_path}"
    try:
        http = get_http()
    except ValueError:
        assert login_required, "Please log in first, use `swanlab login` to log in."
    ds = DataStore()
    ds.open_for_scan(file_path)
    logs = [log for log in ds]
