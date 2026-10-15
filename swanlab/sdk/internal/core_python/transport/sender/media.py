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
from typing import IO, TYPE_CHECKING, Optional, Sequence, Union

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.api.upload import upload_resource
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg.mime import guess_type
from swanlab.sdk.internal.core_python.utils import get_buffer_size
from swanlab.sdk.internal.pkg import adapter, console
from swanlab.sdk.internal.pkg.client.session import SessionWithRetry
from swanlab.sdk.typings.core_python.api.upload import UploadMedia, UploadMediaBatch
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


@dataclass
class MediaUploadPlan:
    """媒体上传收集结果：对象存储载荷、媒体指标与失败记录。"""

    payloads: list[MediaPayload] = field(default_factory=list)
    metrics: UploadMediaBatch = field(default_factory=list)
    metric_records: list[Record] = field(default_factory=list)
    blocked_failures: list[UploadFailure] = field(default_factory=list)


def collect_media_payloads(ctx: CoreContext, records: Sequence[Record]) -> MediaUploadPlan:
    """从 media Record 收集对象存储载荷与媒体指标。"""
    plan = MediaUploadPlan()
    for record in records:
        if not record.HasField("media"):
            plan.blocked_failures.append(
                UploadFailure(source="local_data", message="media record without media payload", records=[record])
            )
            continue
        media_record = record.media
        if not media_record.HasField("timestamp"):
            plan.blocked_failures.append(
                UploadFailure(source="local_data", message="media record without timestamp", records=[record])
            )
            continue
        create_time = media_record.timestamp.ToJsonString()
        metric_chunk: UploadMedia = {
            "key": media_record.key,
            "index": media_record.step,
            "data": [],
            "more": [],
            "create_time": create_time,
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
                plan.payloads.append(
                    MediaPayload(
                        remote_path=remote_path,
                        tracker_key=f"{remote_path}:{size}",
                        display_path=remote_path,
                        size=size,
                        buffer=io.BytesIO(payload),
                        content_type=mime_type,
                    )
                )
                record_paths.append(remote_path)
                if media.caption:
                    metric_chunk["more"].append({"caption": media.caption})
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
            plan.payloads.append(
                MediaPayload(
                    remote_path=remote_path,
                    tracker_key=f"{remote_path}:{size}",
                    display_path=local_path.as_posix(),
                    size=size,
                    buffer=local_path,
                    content_type=mime_type,
                )
            )
            record_paths.append(remote_path)
            if media.caption:
                metric_chunk["more"].append({"caption": media.caption})

        if record_paths:
            metric_chunk["data"] = record_paths
            plan.metrics.append(metric_chunk)
            plan.metric_records.append(record)
        else:
            reason = reasons[0] if reasons else "media record without files"
            plan.blocked_failures.append(UploadFailure(source="local_data", message=reason, records=[record]))
    return plan


def put_media_resources(
    session: SessionWithRetry,
    experiment_id: str,
    payloads: Sequence[MediaPayload],
    *,
    tracker: Optional[UploadTracker] = None,
) -> None:
    """批量上传媒体文件至对象存储，并注册字节进度显示行。"""
    if not payloads:
        return
    for payload in payloads:
        track_file(tracker, payload.tracker_key, payload.display_path, payload.size)
    upload_resource(
        session,
        experiment_id,
        paths=[payload.remote_path for payload in payloads],
        buffers=[payload.buffer for payload in payloads],
        content_types=[payload.content_type for payload in payloads],
        tracker=tracker,
    )


__all__ = [
    "MediaPayload",
    "MediaUploadPlan",
    "collect_media_payloads",
    "put_media_resources",
]
