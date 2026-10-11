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

# ── 媒体载荷收集 ──────────────────────────────────────────────


@dataclass
class MediaPayload:
    """单个媒体文件的对象存储上传载荷。"""

    remote_path: str  # 对象存储路径 media/<type>/<filename>
    tracker_key: str  # tracker 文件显示行注册 key
    display_path: str  # tracker 显示路径（内存载荷用远端路径，本地文件用本地路径）
    size: int
    buffer: Union[IO[bytes], str, Path]  # 内存 buffer 或本地文件路径
    content_type: str


@dataclass
class MediaUploadPlan:
    """一次媒体上传的收集结果：对象存储载荷与媒体指标 DTO。

    payloads 与 metrics、metric_records 一一对应：payload 必然归属某个 metric 的 data 列表，
    metric_records[i] 即产出 metrics[i] 的原 Record；blocked_failures 是未产出任何指标的
    Record 及其原因（缺时间戳、payload 缺失、文件不可用等），由调用方保留并报告。
    """

    payloads: list[MediaPayload] = field(default_factory=list)
    metrics: UploadMediaBatch = field(default_factory=list)
    metric_records: list[Record] = field(default_factory=list)
    blocked_failures: list[UploadFailure] = field(default_factory=list)


def collect_media_payloads(ctx: CoreContext, records: Sequence[Record]) -> MediaUploadPlan:
    """将 media Record 收集为对象存储载荷与媒体指标，纯收集，不发起网络请求。

    - payload 携带内容（skip_store）：内存 buffer 上传，允许空文件，不触碰本地 media 路径
    - skip_store 下缺 payload：告警并归入 blocked，不回退读本地路径
    - 本地文件缺失或读取失败：跳过该文件；整条 Record 均无可用文件时归入 blocked
    """
    plan = MediaUploadPlan()
    for record in records:
        if not record.HasField("media"):
            # 无 media 载荷属确定性脏数据：归入 blocked 保留并报告，不静默消失
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
            # 目前约定的本地文件路径格式为：media/<type>/filename，不区分key
            # 约定保存到对象存储中的文件路径类似于本地文件路径
            medium = adapter.medium[media_record.type]
            remote_path = PurePosixPath("media", medium, media.filename).as_posix()

            if media.HasField("payload"):
                # skip_store：内容在 payload（允许为空文件），直接内存上传，不触碰本地 media 路径
                # media.payload 每次访问都会生成新的 bytes，先取一次用于长度计算和上传流构造
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
                # skip_store 下缺 payload：告警丢弃，不回退读本地路径
                console.warning(f"Media payload missing with skip_store enabled, skipping: {remote_path}")
                reasons.append(f"payload missing with skip_store enabled: {remote_path}")
                continue

            local_path = ctx.media_dir / medium / media.filename
            if not local_path.is_file():
                reasons.append(f"media file not found: {local_path}")
                continue
            try:
                # size 与 mime_type 先于载荷入列完成计算，保证单文件失败时整个 payload
                # 不入列，载荷与指标列表保持原子一致
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


# ── 错误映射 ──────────────────────────────────────────────────


def transport_failure(
    records: Sequence[Record],
    error: BaseException,
    *,
    source: FailureSource,
    stage: Optional[FailureStage] = None,
    message: Optional[str] = None,
) -> UploadFailure:
    """将上传链路异常转为结构化失败：ApiError 提取 HTTP 状态与业务码，其余按异常类型名描述。

    message 缺省时取业务错误信息或异常类型名，均不包含 URL、凭据或完整响应。
    """
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


# ── DTO 编码 ──────────────────────────────────────────────────


def encode_column(record: ColumnRecord) -> UploadColumn:
    """
    将列记录编码为后端所需的格式（DTO）
    """
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


# ── 文件来源解析 ──────────────────────────────────────────────


def resolve_save_source(ctx: CoreContext, save: SaveRecord) -> Optional[Path]:
    """解析 save 记录的可读本地文件路径，原始路径与回退路径均不可读时返回 None。

    ``source_path`` 为绝对路径，跨挂载根 sync 时不可读，回退到 run 目录 ``files`` 子目录：
    CUSTOM 按相对名 ``save.name`` 定位，内部保存按 ``source_path`` 的 basename 定位
    （config 的 ``save.name`` 为 ``"config"`` 而非 ``config.yaml``，故一律用 basename）。
    skip_store 下本 run 不创建镜像，源不可读时直接返回 None（回退命中的是旧 run 遗留文件）。
    """
    primary = Path(save.source_path)
    if primary.is_file():
        return primary
    # skip_store：files 子目录无本 run 镜像，同名文件属于旧 run，回退会上传错误内容
    if ctx.config.skip_store:
        return None
    # 反斜杠路径（在 Windows、sync 在 POSIX）需经 PureWindowsPath 解析分隔符，
    # 本地 Path 才能取到正确的 basename 和子目录层级
    if save.type == SaveType.SAVE_TYPE_CUSTOM:
        fallback = ctx.files_dir / Path(*PureWindowsPath(save.name).parts)
    else:
        fallback = ctx.files_dir / PureWindowsPath(save.source_path).name
    if fallback.is_file():
        console.debug(f"Save source path not readable, recovered from run directory: {fallback}")
        return fallback
    return None


# ── hash 与文件标识 ───────────────────────────────────────────


def compute_md5(source: Union[str, bytes], chunk_size: int = 8 * 1024 * 1024) -> str:
    """计算 MD5。接受文件路径（str）或已读取的字节数据（bytes）。"""
    md5 = hashlib.md5()
    if isinstance(source, bytes):
        md5.update(source)
    else:
        with open(source, "rb") as f:
            for chunk in iter(lambda: f.read(chunk_size), b""):
                md5.update(chunk)
    return md5.hexdigest()


def save_tracker_key(file_info: SaveFileEntry) -> str:
    """获取 save 文件在 tracker 字典中的注册 key，以远端 path 作为稳定前缀。"""
    return f"{file_info['path']}:{file_info['size']}:{file_info['md5']}"


__all__ = [
    "MediaPayload",
    "MediaUploadPlan",
    "collect_media_payloads",
    "compute_md5",
    "encode_column",
    "resolve_save_source",
    "save_tracker_key",
    "transport_failure",
]
