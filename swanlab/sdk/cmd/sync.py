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
    SyncRange,
)
from swanlab.proto.swanlab.settings.core.v1.core_pb2 import CoreSettings
from swanlab.sdk.cmd import utils
from swanlab.sdk.cmd.login import login_raw
from swanlab.sdk.internal import impl
from swanlab.sdk.internal.core_python import client
from swanlab.sdk.internal.pkg import console, helper
from swanlab.sdk.internal.run.progress import run_with_progress
from swanlab.sdk.internal.settings import Settings, create_settings
from swanlab.sdk.protocol.core import CoreSyncProtocol
from swanlab.utils.experiment import generate_id
from swanlab.utils.time import TimeLike, format_time, set_timestamp

__all__ = ["sync", "ensure_run_dir"]


def sync(
    run_dir: Union[Path, str],
    settings: Optional[Settings] = None,
    *,
    step_start: Optional[int] = None,
    step_end: Optional[int] = None,
    start_time: Optional[TimeLike] = None,
    end_time: Optional[TimeLike] = None,
):
    """Synchronize local SwanLab run data to the cloud.

    Uploads a local SwanLab run directory, usually created in offline mode,
    to SwanLab Cloud. Authentication must be available before syncing, either
    from existing credentials or from the API key configured in `settings`.

    Passing any of ``step_start``/``step_end``/``start_time``/``end_time`` enables
    range sync: only scalar/media records inside the step and time window, plus
    terminal logs inside the time window, are replayed. Range sync bypasses the
    duplicate-step guard, so it can backfill history into an existing experiment
    by passing ``settings.run.id``.

    Sync always creates or resumes the target experiment first. A finished target
    is resumed to RUNNING, so the replayed rows stay inside the server-side query
    window; afterwards the local FinishRecord decides the final state (a new
    experiment without a FinishRecord is reported as CRASHED, while an existing
    experiment without one is left RUNNING for other writers). Range bounds are
    inclusive.

    .. note::
        A step-only range does not replay terminal logs; pass ``start_time`` /
        ``end_time`` to include them.

    :param run_dir: Existing local SwanLab run directory to synchronize.
    :param settings: Optional sync configuration, including API key, host, workspace, project, and run options. Uses default Settings if not provided.
    :param step_start: Only sync scalar/media records whose step is >= this value.
    :param step_end: Only sync scalar/media records whose step is <= this value.
    :param start_time: Only sync records at or after this time. Times without an offset
        are interpreted in the system local timezone.
    :param end_time: Only sync records at or before this time. Same as ``start_time``.
    """
    if not isinstance(run_dir, Path):
        run_dir = Path(run_dir)
    # 0. 构建范围过滤器，参数非法时直接失败，避免创建实验后才发现问题
    sync_range = build_sync_range(
        step_start=step_start,
        step_end=step_end,
        start_time=start_time,
        end_time=end_time,
    )
    # step 范围只约束指标，终端日志仅受时间范围约束
    if sync_range is not None and not (sync_range.HasField("start_time") or sync_range.HasField("end_time")):
        console.warning(
            "Terminal logs are not replayed for a step-only range; pass start_time/end_time to include them."
        )
    # 1. settings 合并，并做业务层额外参数限制
    sync_settings = create_settings()
    if settings is not None:
        sync_settings.merge_settings(settings)
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
    assert client.exists()

    if sync_range is not None:
        console.info(f"Syncing records in {format_sync_range(sync_range)}")

    # 3. 初始化core sync对象，并启动sync
    core = impl.create_core_sync()
    flush_resp = _deliver_sync_start(core, core_settings, workspace, project, run_id, sync_range)

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


def build_sync_range(
    *,
    step_start: Optional[int] = None,
    step_end: Optional[int] = None,
    start_time: Optional[TimeLike] = None,
    end_time: Optional[TimeLike] = None,
) -> Optional[SyncRange]:
    """构建范围同步协议对象，所有参数为空时返回 None，表示全量同步。

    不带时区的时间按系统本地时区解释。
    """
    if step_start is None and step_end is None and start_time is None and end_time is None:
        return None
    if step_start is not None and step_start < 0:
        raise ValueError("step_start must be a non-negative integer")
    if step_end is not None and step_end < 0:
        raise ValueError("step_end must be a non-negative integer")
    if step_start is not None and step_end is not None and step_start > step_end:
        raise ValueError("step_start must be less than or equal to step_end")

    sync_range = SyncRange()
    if step_start is not None:
        sync_range.start_step = step_start
    if step_end is not None:
        sync_range.end_step = step_end
    if start_time is not None:
        set_timestamp(sync_range.start_time, start_time)
    if end_time is not None:
        set_timestamp(sync_range.end_time, end_time)
    if start_time is not None and end_time is not None:
        if sync_range.start_time.ToNanoseconds() > sync_range.end_time.ToNanoseconds():
            raise ValueError("start_time must be earlier than or equal to end_time")
    return sync_range


def format_sync_range(sync_range: SyncRange) -> str:
    """将范围同步配置格式化为可读文本，用于终端提示（按系统本地时区展示）。"""
    parts = []
    if sync_range.HasField("start_step") or sync_range.HasField("end_step"):
        start = str(sync_range.start_step) if sync_range.HasField("start_step") else "*"
        end = str(sync_range.end_step) if sync_range.HasField("end_step") else "*"
        parts.append(f"step {start}..{end}")
    if sync_range.HasField("start_time") or sync_range.HasField("end_time"):
        start = format_time(sync_range.start_time) if sync_range.HasField("start_time") else "*"
        end = format_time(sync_range.end_time) if sync_range.HasField("end_time") else "*"
        parts.append(f"time {start}..{end}")
    return ", ".join(parts)


@utils.with_loading_animation("Starting core to sync data...")
def _deliver_sync_start(
    core: CoreSyncProtocol,
    core_settings: CoreSettings,
    workspace: Optional[str],
    project: Optional[str],
    run_id: Optional[str],
    sync_range: Optional[SyncRange] = None,
) -> DeliverSyncFlushResponse:
    start_response = core.deliver_sync_start(
        DeliverSyncStartRequest(
            core_settings=core_settings,
            workspace=workspace,
            project=project,
            id=run_id,
            range=sync_range,
        )
    )
    if not start_response.success:
        raise RuntimeError(start_response.message)
    resp = core.deliver_sync_flush()
    if not resp.success:
        raise RuntimeError(resp.message)
    return resp
