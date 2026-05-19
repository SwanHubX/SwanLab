"""
@author: cunyue
@file: sync.py
@time: 2026/5/13 15:56
@description: 同步实验数据到云端
"""

import os
from pathlib import Path
from typing import Optional, Union

from pydantic import DirectoryPath, TypeAdapter
from rich.text import Text

from swanlab.exceptions import AuthenticationError
from swanlab.proto.swanlab.grpc.core.v1.sync_pb2 import (
    DeliverSyncFlushResponse,
    DeliverSyncStartRequest,
)
from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings
from swanlab.sdk.cmd import utils
from swanlab.sdk.cmd.login import login_raw
from swanlab.sdk.internal import impl
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import console, helper
from swanlab.sdk.internal.run.progress import run_with_progress
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings import settings as global_settings
from swanlab.sdk.protocol.core import CoreSyncProtocol
from swanlab.utils.experiment import generate_id

__all__ = ["sync", "ensure_run_dir"]


def sync(run_dir: Union[Path, str], settings: Optional[Settings] = None):
    """Synchronize local SwanLab run data to the cloud.

    Uploads a local SwanLab run directory, usually created in offline mode,
    to SwanLab Cloud. Authentication must be available before syncing, either
    from existing credentials or from the API key configured in `settings`.

    :param run_dir: Existing local SwanLab run directory to synchronize.
    :param settings: Optional sync configuration, including API key, host, workspace, project, and run options. Uses default Settings if not provided.
    """
    if not isinstance(run_dir, Path):
        run_dir = Path(run_dir)
    # 1. settings 合并，并做业务层额外参数限制
    sync_settings = Settings()
    if settings is not None:
        sync_settings.merge_settings(settings)
    sync_settings.merge_settings(global_settings)
    workspace = sync_settings.project.workspace
    project = sync_settings.project.name
    run_dir = ensure_run_dir(run_dir)
    run_id = sync_settings.run.id or generate_id()
    core_settings = sync_settings.to_core_proto(run_id, run_dir)

    # 2. 校验 api key，启动client
    if not client.exists():
        if sync_settings.api_key is None:
            raise AuthenticationError("Please login first, or use `swanlab.login()`")
        login_raw(sync_settings.api_key, host=sync_settings.api_host, print_welcome=False)
        sync_settings.merge_settings(global_settings)
    assert client.exists()

    # 3. 初始化core sync对象，并启动sync
    core = impl.create_core_sync()
    flush_resp = _deliver_sync_start(core, core_settings, workspace, project, run_id)

    # 4. 打印创建的实验 URL
    if flush_resp and flush_resp.path:
        run_path = helper.fmt_run_path(flush_resp.path)
        web_host = sync_settings.web_host
        run_url = f"{web_host}{run_path}"
        project_url = run_url.split("/runs/")[0]
        console.info("📁 View project at", Text(project_url, style=f"link {project_url} blue underline"))
        console.info("🚀 View synced run at", Text(run_url, style=f"link {run_url} blue underline"))

    # 5. 等待core完成同步
    confirm_resp = run_with_progress(
        stats_fn=core.get_operation_stats,
        blocking_fn=core.confirm_sync_finish,
        unit="auto",
    )

    if not confirm_resp.success:
        console.error(f"Confirm sync finish failed: {confirm_resp.message}")


def ensure_run_dir(run_dir: DirectoryPath) -> DirectoryPath:
    """校验 run_dir 存在、是目录、具备可读权限，并返回绝对路径。"""
    run_dir = TypeAdapter(DirectoryPath).validate_python(run_dir)
    run_dir = run_dir.resolve()
    if not os.access(run_dir, os.R_OK):
        raise PermissionError(f"SwanLab run directory is not readable: {run_dir}")
    return run_dir


@utils.with_loading_animation("Starting core to sync data...")
def _deliver_sync_start(
    core: CoreSyncProtocol,
    core_settings: CoreSettings,
    workspace: Optional[str],
    project: Optional[str],
    run_id: Optional[str],
) -> DeliverSyncFlushResponse:
    start_response = core.deliver_sync_start(
        DeliverSyncStartRequest(
            core_settings=core_settings,
            workspace=workspace,
            project=project,
            id=run_id,
        )
    )
    if not start_response.success:
        raise RuntimeError(start_response.message)
    resp = core.deliver_sync_flush()
    if not resp.success:
        raise RuntimeError(resp.message)
    return resp
