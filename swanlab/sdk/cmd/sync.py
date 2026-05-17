"""
@author: cunyue
@file: sync.py
@time: 2026/5/13 15:56
@description: 同步实验数据到云端
"""

from pathlib import Path
from typing import Optional

from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.protocol.core import CoreSyncProtocol


def sync(run_dir: Path, settings: Optional[Settings] = None):
    """Synchronize local SwanLab run data to the cloud.

    Uploads a local SwanLab run directory, usually created in offline mode,
    to SwanLab Cloud. Authentication must be available before syncing, either
    from existing credentials or from the API key configured in `settings`.

    :param run_dir: Local SwanLab run directory to synchronize.
    :param settings: Optional sync configuration, including API key, host, workspace, project, and run options. Uses default Settings if not provided.
    """
    # 1. settings 合并，并做业务层额外参数限制
    sync_settings = Settings()
    if settings is not None:
        sync_settings.merge_settings(settings)
    # project 参数必须存在
    if sync_settings.project.name is None:
        raise ValueError("project.name is required when use sync")
    # 2. 初始化core sync对象
    _ = _create_core_sync()


def _create_core_sync() -> CoreSyncProtocol:
    return CoreSyncPython()
