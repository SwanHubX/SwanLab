"""
@author: caddiesnew
@file: media.py
@time: 2026/4/16
@description: 媒体载荷收集与对象存储传输
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import IO, TYPE_CHECKING, Dict, Optional, Sequence, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.api.upload import upload_resource
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg.mime import guess_type
from swanlab.sdk.internal.core_python.utils import get_buffer_size
from swanlab.sdk.internal.pkg import adapter, console
from swanlab.sdk.internal.pkg.client.session import SessionWithRetry
from swanlab.sdk.typings.core_python.api.upload import UploadMedia
from swanlab.sdk.typings.core_python.transport import UploadFailure

from .utils import track_file

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


@dataclass
class MediaPayload:
    """单个媒体文件的上传载荷。"""

    remote_path: str
    tracker_key: str
    display_path: str
    size: int
    buffer: Union[IO[bytes], str, Path]
    content_type: str
    # 归属上传任务的下标，对应 MediaRecordTask.index
    task_index: int


@dataclass
class MediaRecordTask:
    """一条 media record 的上传任务：媒体指标与逐文件载荷。"""

    # 扁平载荷列表中的归属下标，对应 MediaPayload.task_index
    index: int
    record: Record
    payloads: list[MediaPayload] = field(default_factory=list)
    # 无可传文件时为 None，此时 failure 记录终结原因
    metric: Optional[UploadMedia] = None
    # 与 metric 互斥：已终结的任务携带失败诊断
    failure: Optional[UploadFailure] = None


def collect_media_tasks(ctx: CoreContext, records: Sequence[Record]) -> list[MediaRecordTask]:
    """把 media Record 拆成逐记录的上传任务：指标、载荷与终结原因。"""
    tasks: list[MediaRecordTask] = []
    for record in records:
        task = MediaRecordTask(index=len(tasks), record=record)
        tasks.append(task)
        if not record.HasField("media"):
            task.failure = UploadFailure(message="media record without media payload", records=[record])
            continue
        media_record = record.media
        if not media_record.HasField("timestamp"):
            task.failure = UploadFailure(message="media record without timestamp", records=[record])
            continue
        metric: UploadMedia = {
            "key": media_record.key,
            "index": media_record.step,
            "data": [],
            "more": [],
            "create_time": media_record.timestamp.ToJsonString(),
        }
        record_paths: list[str] = []
        reasons: list[str] = []
        for media in media_record.value.items:
            # 路径格式：media/<type>/<filename>
            medium = adapter.medium[media_record.type]
            remote_path = PurePosixPath("media", medium, media.filename).as_posix()

            if media.HasField("payload"):
                payload = media.payload
                size = len(payload)
                mime_type = guess_type(media.filename)
                item_payload = MediaPayload(
                    remote_path=remote_path,
                    tracker_key=f"{remote_path}:{size}",
                    display_path=remote_path,
                    size=size,
                    buffer=io.BytesIO(payload),
                    content_type=mime_type,
                    task_index=task.index,
                )
                task.payloads.append(item_payload)
                record_paths.append(remote_path)
                if media.caption:
                    metric["more"].append({"caption": media.caption})
                continue

            if ctx.config.skip_store:
                console.warning(f"Media payload missing with skip_store enabled, skipping: {remote_path}")
                reasons.append(f"payload missing with skip_store enabled: {remote_path}")
                continue

            local_path = ctx.media_dir / medium / media.filename
            if not local_path.is_file():
                reasons.append(f"media file not found: {local_path}")
                continue
            try:
                size = get_buffer_size(local_path)
                mime_type = guess_type(local_path)
            except OSError as e:
                reasons.append(f"media file unreadable: {local_path}: {e}")
                continue
            item_payload = MediaPayload(
                remote_path=remote_path,
                tracker_key=f"{remote_path}:{size}",
                display_path=local_path.as_posix(),
                size=size,
                buffer=local_path,
                content_type=mime_type,
                task_index=task.index,
            )
            task.payloads.append(item_payload)
            record_paths.append(remote_path)
            if media.caption:
                metric["more"].append({"caption": media.caption})

        if record_paths:
            metric["data"] = record_paths
            task.metric = metric
        else:
            reason = reasons[0] if reasons else "media record without files"
            task.failure = UploadFailure(message=reason, records=[record])
    return tasks


def put_media_resources(
    session: SessionWithRetry,
    experiment_id: str,
    payloads: Sequence[MediaPayload],
    *,
    tracker: Optional[UploadTracker] = None,
) -> Dict[int, BaseException]:
    """批量上传媒体文件至对象存储，并注册字节进度显示行。返回 {索引: 异常} 字典。"""
    if not payloads:
        return {}
    for payload in payloads:
        track_file(tracker, payload.tracker_key, payload.display_path, payload.size)
    return upload_resource(
        session,
        experiment_id,
        paths=[payload.remote_path for payload in payloads],
        buffers=[payload.buffer for payload in payloads],
        content_types=[payload.content_type for payload in payloads],
        tracker=tracker,
    )


__all__ = [
    "MediaPayload",
    "MediaRecordTask",
    "collect_media_tasks",
    "put_media_resources",
]
