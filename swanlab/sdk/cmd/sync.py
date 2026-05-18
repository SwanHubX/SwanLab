"""
@author: cunyue
@file: sync.py
@time: 2026/5/13 15:56
@description: 同步实验数据到云端
"""

import os
from typing import Optional

from pydantic import DirectoryPath, TypeAdapter

from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import DeliverSyncStartRequest
from swanlab.sdk.internal.core_python.sync import CoreSyncPython
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.protocol.core import CoreSyncProtocol
from swanlab.utils.experiment import generate_id

__all__ = ["sync", "ensure_run_dir"]


def sync(run_dir: DirectoryPath, settings: Optional[Settings] = None):
    """Synchronize local SwanLab run data to the cloud.

    Uploads a local SwanLab run directory, usually created in offline mode,
    to SwanLab Cloud. Authentication must be available before syncing, either
    from existing credentials or from the API key configured in `settings`.

    :param run_dir: Existing local SwanLab run directory to synchronize.
    :param settings: Optional sync configuration, including API key, host, workspace, project, and run options. Uses default Settings if not provided.
    """
    # 1. settings 合并，并做业务层额外参数限制
    sync_settings = Settings()
    if settings is not None:
        sync_settings.merge_settings(settings)
    project = sync_settings.project.name
    run_dir = ensure_run_dir(run_dir)
    run_id = sync_settings.run.id or generate_id()
    core_settings = sync_settings.to_core_proto(run_id, run_dir)

    # 2. 初始化core sync对象，并启动sync
    core = _create_core_sync()
    start_response = core.deliver_sync_start(
        DeliverSyncStartRequest(
            core_settings=core_settings,
            project=project,
            id=run_id,
        )
    )
    if not start_response.success:
        raise RuntimeError(start_response.message)
    # TODO: 3. 等待core完全读取所有的record

    # 4. 告诉core sync同步结束，并且sdk等待core完成
    core.deliver_sync_finish()


def ensure_run_dir(run_dir: DirectoryPath) -> DirectoryPath:
    """校验 run_dir 存在、是目录、具备可读权限，并返回绝对路径。"""
    run_dir = TypeAdapter(DirectoryPath).validate_python(run_dir)
    run_dir = run_dir.resolve()
    if not os.access(run_dir, os.R_OK):
        raise PermissionError(f"SwanLab run directory is not readable: {run_dir}")
    return run_dir


def _create_core_sync() -> CoreSyncProtocol:
    if helper.get_core_impl() == "python":
        return CoreSyncPython()
    else:
        # TODO: Core 微服务无感接入
        raise NotImplementedError("The SwanLab Go core sync runtime is not available yet.")
