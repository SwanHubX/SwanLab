#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@author: caddiesnew
@file: utils.py
@time: 2026/3/30 15:33
@description: sender 辅助函数：DTO 编码、媒体载荷收集、文件来源解析、hash 与文件标识
"""

import hashlib
import io
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import IO, Optional, Sequence, Union

import requests

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg.mime import guess_type
from swanlab.sdk.internal.core_python.utils import get_buffer_size
from swanlab.sdk.internal.pkg import adapter, console
from swanlab.sdk.typings.core_python.api.save import SaveFileEntry
from swanlab.sdk.typings.core_python.api.upload import UploadColumn, UploadMedia, UploadMediaBatch
from swanlab.sdk.typings.core_python.transport import FailureSource, FailureStage, UploadFailure


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


def transport_failure(
    records: Sequence[Record],
    error: BaseException,
    *,
    source: FailureSource,
    stage: Optional[FailureStage] = None,
    message: Optional[str] = None,
) -> UploadFailure:
    """将上传异常转换为 UploadFailure。"""
    if isinstance(error, ApiError):
        return UploadFailure(
            source=source,
            stage=stage,
            message=message if message is not None else str(error.message),
            records=list(records),
            status=getattr(error.response, "status_code", None),
            code=error.code,
            exception=error,
        )
    return UploadFailure(
        source=source,
        stage=stage,
        message=message if message is not None else type(error).__name__,
        records=list(records),
        exception=error,
    )


def is_retryable_error(error: requests.RequestException) -> bool:
    """判断请求失败是否属于可恢复的网络或服务端故障。"""
    if not isinstance(error, ApiError):
        return True
    status = getattr(error.response, "status_code", None)
    return status is None or status >= 500


def encode_column(record: ColumnRecord) -> UploadColumn:
    """将列记录编码为后端所需的格式（DTO）。"""
    # TODO 在 v0.11 版本 精简一下 column 的类型定义
    column: UploadColumn = {"key": record.column_key, "type": adapter.column[record.column_type]}
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["type"] = "SYSTEM"
    if record.HasField("x_axis"):
        column["xAxis"] = record.x_axis
    if record.section_name:
        column["sectionName"] = record.section_name
    if record.hidden:
        column["hidden"] = True
    return column


def resolve_save_source(ctx: CoreContext, save: SaveRecord) -> Optional[Path]:
    """解析 save 记录的本地文件路径，源路径不可读时回退至 run 目录下的 files 子目录。"""
    primary = Path(save.source_path)
    if primary.is_file():
        return primary
    if ctx.config.skip_store:
        return None
    if save.type == SaveType.SAVE_TYPE_CUSTOM:
        fallback = ctx.files_dir / Path(*PureWindowsPath(save.name).parts)
    else:
        fallback = ctx.files_dir / PureWindowsPath(save.source_path).name
    if fallback.is_file():
        console.debug(f"Save source path not readable, recovered from run directory: {fallback}")
        return fallback
    return None


def compute_md5(source: Union[str, bytes], chunk_size: int = 8 * 1024 * 1024) -> str:
    """计算文件路径或字节数据的 MD5。"""
    md5 = hashlib.md5()
    if isinstance(source, bytes):
        md5.update(source)
    else:
        with open(source, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                md5.update(chunk)
    return md5.hexdigest()


def save_tracker_key(file_info: SaveFileEntry) -> str:
    """构造 save 文件在 tracker 中的注册 key。"""
    return f"{file_info['path']}:{file_info['size']}:{file_info['md5']}"


__all__ = [
    "MediaPayload",
    "MediaUploadPlan",
    "collect_media_payloads",
    "compute_md5",
    "encode_column",
    "is_retryable_error",
    "resolve_save_source",
    "save_tracker_key",
    "transport_failure",
]
