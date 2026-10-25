"""
@author: cunyue
@file: sync.py
@time: 2026/5/13 15:56
@description: 同步实验数据到云端
"""

import os
from datetime import datetime, timezone, tzinfo
from pathlib import Path
from typing import Optional, Union

from google.protobuf.timestamp_pb2 import Timestamp
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
from swanlab.utils.time import parse_timestamp_ns, parse_timezone

__all__ = ["sync", "ensure_run_dir"]

TimeLike = Union[int, float, str, datetime]
"""sync 范围同步允许的时间输入：Unix 时间戳（秒/毫秒）、ISO 8601 字符串或 datetime。"""

TimezoneLike = Union[str, tzinfo]
"""sync 范围同步允许的时区输入：IANA 时区名、UTC 偏移或 tzinfo。"""


def sync(
    run_dir: Union[Path, str],
    settings: Optional[Settings] = None,
    *,
    step_start: Optional[int] = None,
    step_end: Optional[int] = None,
    start_time: Optional[TimeLike] = None,
    end_time: Optional[TimeLike] = None,
    timezone: Optional[TimezoneLike] = None,
):
    """Synchronize local SwanLab run data to the cloud.

    Uploads a local SwanLab run directory, usually created in offline mode,
    to SwanLab Cloud. Authentication must be available before syncing, either
    from existing credentials or from the API key configured in `settings`.

    Passing any of ``step_start``/``step_end``/``start_time``/``end_time`` enables
    range sync: only scalar/media records inside the step and time window, plus
    terminal logs inside the time window, are replayed. Range sync bypasses the
    duplicate-step guard, so it can backfill history into an existing experiment
    by passing ``settings.run.id`` (in-place restore); in that case the experiment
    lifecycle state is left untouched. Range bounds are inclusive.

    :param run_dir: Existing local SwanLab run directory to synchronize.
    :param settings: Optional sync configuration, including API key, host, workspace, project, and run options. Uses default Settings if not provided.
    :param step_start: Only sync scalar/media records whose step is >= this value.
    :param step_end: Only sync scalar/media records whose step is <= this value.
    :param start_time: Only sync records at or after this time.
    :param end_time: Only sync records at or before this time.
    :param timezone: Timezone used to interpret and display times without an explicit
        offset, e.g. ``"Asia/Shanghai"`` or ``"+08:00"``. Defaults to the system timezone.
        Inputs that already carry an offset (e.g. ``"2026-10-20T10:00:00Z"``) keep their own.
    """
    if not isinstance(run_dir, Path):
        run_dir = Path(run_dir)
    # 0. 构建范围过滤器，参数非法时直接失败，避免创建实验后才发现问题
    tz = parse_timezone(timezone)
    sync_range = build_sync_range(
        step_start=step_start,
        step_end=step_end,
        start_time=start_time,
        end_time=end_time,
        timezone=tz,
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
        console.info(f"Syncing records in {format_sync_range(sync_range, timezone=tz)}")

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
    timezone: Optional[TimezoneLike] = None,
) -> Optional[SyncRange]:
    """构建范围同步协议对象，所有参数为空时返回 None，表示全量同步。"""
    if step_start is None and step_end is None and start_time is None and end_time is None:
        return None
    if step_start is not None and step_start < 0:
        raise ValueError("step_start must be a non-negative integer")
    if step_end is not None and step_end < 0:
        raise ValueError("step_end must be a non-negative integer")
    if step_start is not None and step_end is not None and step_start > step_end:
        raise ValueError("step_start must be less than or equal to step_end")

    tz = parse_timezone(timezone)
    sync_range = SyncRange()
    if step_start is not None:
        sync_range.start_step = step_start
    if step_end is not None:
        sync_range.end_step = step_end
    if start_time is not None:
        _set_timestamp(sync_range.start_time, start_time, tz)
    if end_time is not None:
        _set_timestamp(sync_range.end_time, end_time, tz)
    if start_time is not None and end_time is not None:
        if sync_range.start_time.ToNanoseconds() > sync_range.end_time.ToNanoseconds():
            raise ValueError("start_time must be earlier than or equal to end_time")
    return sync_range


def format_sync_range(sync_range: SyncRange, *, timezone: Optional[TimezoneLike] = None) -> str:
    """将范围同步配置格式化为可读文本，用于终端提示。"""
    tz = parse_timezone(timezone)
    parts = []
    if sync_range.HasField("start_step") or sync_range.HasField("end_step"):
        start = str(sync_range.start_step) if sync_range.HasField("start_step") else "*"
        end = str(sync_range.end_step) if sync_range.HasField("end_step") else "*"
        parts.append(f"step {start}..{end}")
    if sync_range.HasField("start_time") or sync_range.HasField("end_time"):
        start = _format_time(sync_range.start_time, tz) if sync_range.HasField("start_time") else "*"
        end = _format_time(sync_range.end_time, tz) if sync_range.HasField("end_time") else "*"
        parts.append(f"time {start}..{end}")
    return ", ".join(parts)


def _set_timestamp(target: Timestamp, value: TimeLike, tz: Optional[tzinfo]) -> None:
    try:
        target.FromNanoseconds(parse_timestamp_ns(value, tz=tz))
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid time value {value!r}: {e}") from e


def _format_time(ts: Timestamp, tz: Optional[tzinfo]) -> str:
    dt = ts.ToDatetime().replace(tzinfo=timezone.utc)
    return (dt.astimezone(tz) if tz is not None else dt.astimezone()).strftime("%Y-%m-%d %H:%M:%S")


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
