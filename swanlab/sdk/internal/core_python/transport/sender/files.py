"""
@author: caddiesnew
@file: files.py
@time: 2026/4/16
@description: 媒体与 CUSTOM save 的对象存储传输：presigned PUT 并发、分片上传、字节进度上报

本模块是自由函数集合，session/tracker 等依赖全部显式传参，不持有模块级状态；
Record 编码与收集见 utils.py，record 路由与业务指标提交见 sender/__init__.py。
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal, Optional, Protocol, cast

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.api.save import (
    complete_multipart_save,
    complete_save_files,
    prepare_multipart_save,
    prepare_save_files,
)
from swanlab.sdk.internal.core_python.api.upload import upload_resource, upload_saves
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg.mime import guess_type
from swanlab.sdk.internal.core_python.utils import ProgressFileWrapper
from swanlab.sdk.internal.pkg import console, safe
from swanlab.sdk.internal.pkg.client.session import SessionWithRetry
from swanlab.sdk.internal.pkg.executor import SafeThreadPoolExecutor
from swanlab.sdk.typings.core_python.api.save import (
    CompletedMultipartSaveFile,
    CompletedPart,
    CompleteSaveFile,
    MultipartPart,
    PrepareMultipartSaveFile,
    PrepareSaveFile,
    SaveFileEntry,
)
from swanlab.sdk.typings.core_python.api.upload import UploadResource

from .utils import MediaPayload, compute_md5, resolve_save_source, save_tracker_key

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


class PartProgressReporter(Protocol):
    """分片字节进度上报回调：(文件 key, 分片 upload key, 当前字节数, final 是否完成)。"""

    def __call__(self, file_key: str, upload_key: str, current_bytes: int, *, final: bool = False) -> None: ...


def _track_file(tracker: Optional[UploadTracker], key: str, path: str, size: int) -> None:
    """在 tracker 中注册文件显示行，重复 key 由 UploadTracker 自动去重。"""
    if tracker is None or size <= 0:
        return
    tracker.add_file(key=key, path=path, total=size)


# ── 媒体 ──────────────────────────────────────────────────────


def put_media_resources(
    session: SessionWithRetry,
    experiment_id: str,
    payloads: Sequence[MediaPayload],
    *,
    tracker: Optional[UploadTracker] = None,
) -> None:
    """媒体载荷批量传输到对象存储（presigned PUT 并发），并注册字节进度显示行。"""
    if not payloads:
        return
    for payload in payloads:
        _track_file(tracker, payload.tracker_key, payload.display_path, payload.size)
    upload_resource(
        session,
        experiment_id,
        paths=[payload.remote_path for payload in payloads],
        buffers=[payload.buffer for payload in payloads],
        content_types=[payload.content_type for payload in payloads],
        tracker=tracker,
    )


# ── CUSTOM save ───────────────────────────────────────────────


def upload_save_files(
    ctx: CoreContext,
    session: SessionWithRetry,
    records: Sequence[Record],
    *,
    tracker: Optional[UploadTracker] = None,
    on_part_progress: Optional[PartProgressReporter] = None,
) -> None:
    """CUSTOM save 文件上传编排：解析来源 → 校验大小 → 按 save_batch 分批 → 小/大文件分流。

    调用方保证 records 均为 CUSTOM 类型且不携带 payload。
    """
    config = ctx.config
    # 1. 收集合法文件（不读内容，避免大批次内存爆炸）
    pending: list[tuple[str, int, str]] = []  # (source_path, size, name)
    for record in records:
        if not record.HasField("save"):
            continue
        save = record.save
        source_ref_path = resolve_save_source(ctx, save)
        if source_ref_path is None:
            console.warning(f"Save file not found, skipping: {save.source_path}")
            continue
        size: Optional[int] = None
        with safe.block(message=f"Failed to stat save file, skipping: {source_ref_path}"):
            size = source_ref_path.stat().st_size
        if size is None:
            continue
        if size > config.save_size:
            console.warning(f"Save file exceeds size limit ({size} > {config.save_size}), skipping: {source_ref_path}")
            continue
        # 用解析后的可读路径（可能是回退到 files 子目录的路径）作为实际读取路径
        # 保持 str() 原生路径形式，与下游 compute_md5 / upload_saves 返回的 source_path 字符串一致，
        # 避免 Windows 上 .as_posix() 与 str() 形式不一致导致 source_to_path 映射查不到
        pending.append((str(source_ref_path), size, save.name))
    if not pending:
        console.warning("No valid files to save.")
        return
    if config.save_batch <= 0:
        console.warning(f"Invalid save batch size ({config.save_batch}), skipping save upload")
        return

    # 2. 按 save_batch 拆分文件列表，避免单次 prepare/complete 超过后端限制
    for index in range(0, len(pending), config.save_batch):
        _upload_save_batch(
            ctx,
            session,
            pending[index : index + config.save_batch],
            tracker=tracker,
            on_part_progress=on_part_progress,
        )


def _upload_save_batch(
    ctx: CoreContext,
    session: SessionWithRetry,
    pending: Sequence[tuple[str, int, str]],
    *,
    tracker: Optional[UploadTracker] = None,
    on_part_progress: Optional[PartProgressReporter] = None,
) -> None:
    """上传一个 save_batch 内的文件列表。"""
    config = ctx.config
    with safe.block(message="Failed to upload save files, skipping"):
        file_entries: list[SaveFileEntry] = []
        with SafeThreadPoolExecutor(max_workers=4) as executor:
            md5_futures = {source_path: executor.submit(compute_md5, source_path) for source_path, _, _ in pending}
            for source_path, size, name in pending:
                md5 = md5_futures[source_path].result()
                if not isinstance(md5, str):
                    console.warning(f"Save file MD5 computation failed, skipping: {source_path}")
                    continue
                mime_type = guess_type(source_path)
                file_entries.append(
                    {
                        "path": name,
                        "source_path": source_path,
                        # make pycharm happy
                        "size": cast(int, size),
                        "md5": md5,
                        "mime_type": mime_type,
                    }
                )

        small_files = [f for f in file_entries if f["size"] < config.save_split]
        large_files = [f for f in file_entries if f["size"] >= config.save_split]

        if small_files:
            _upload_small_saves(ctx, session, small_files, tracker=tracker)
        if large_files:
            _upload_large_saves(ctx, session, large_files, tracker=tracker, on_part_progress=on_part_progress)


def _upload_small_saves(
    ctx: CoreContext,
    session: SessionWithRetry,
    files: list[SaveFileEntry],
    *,
    tracker: Optional[UploadTracker] = None,
) -> None:
    """小文件：prepare → presigned PUT（并发）→ complete（逐文件状态）。"""
    prepare_items: list[PrepareSaveFile] = [
        {"path": f["path"], "size": f["size"], "md5": f["md5"], "mimeType": f["mime_type"]} for f in files
    ]
    resp = prepare_save_files(ctx.experiment_id, files=prepare_items)
    urls = resp["urls"]
    if len(urls) != len(files):
        console.warning(f"Save prepare returned {len(urls)} URLs for {len(files)} files, skipping")
        return
    resources: list[UploadResource] = []
    for f, url in zip(files, urls):
        tracker_key = save_tracker_key(f)
        _track_file(tracker, tracker_key, f["source_path"], f["size"])
        resources.append(
            {
                "url": url,
                "source_path": f["source_path"],
                "content_type": f["mime_type"],
                "size": f["size"],
                "tracker_key": tracker_key,
            }
        )
    source_to_path = {f["source_path"]: f["path"] for f in files}
    statuses: dict[str, str] = {f["path"]: "FAILED" for f in files}
    with safe.block(message="Failed to upload save files, skipping"):
        uploaded_paths = upload_saves(
            session,
            resources=resources,
            tracker=tracker,
        )
        for source_path in uploaded_paths:
            statuses[source_to_path[source_path]] = "UPLOADED"

    complete_items: list[CompleteSaveFile] = [
        {"path": path, "state": cast(Literal["UPLOADED", "FAILED"], state)} for path, state in statuses.items()
    ]
    with safe.block(message="Failed to complete save upload"):
        complete_save_files(ctx.experiment_id, files=complete_items)
    synced = sum(1 for s in statuses.values() if s == "UPLOADED")
    if synced > 0:
        console.info(f"Synced {synced} save file(s)")


def _upload_large_saves(
    ctx: CoreContext,
    session: SessionWithRetry,
    files: list[SaveFileEntry],
    *,
    tracker: Optional[UploadTracker] = None,
    on_part_progress: Optional[PartProgressReporter] = None,
) -> None:
    """大文件：分片上传。"""
    part_size = ctx.config.save_part
    prepare_items: list[PrepareMultipartSaveFile] = []
    for f in files:
        part_count = math.ceil(f["size"] / part_size)
        prepare_items.append(
            {
                "path": f["path"],
                "size": f["size"],
                "md5": f["md5"],
                "mimeType": f["mime_type"],
                "count": part_count,
            }
        )
    resp = prepare_multipart_save(ctx.experiment_id, files=prepare_items)
    multipart_files = resp["files"]
    if len(multipart_files) != len(files):
        console.warning(f"Multipart prepare returned {len(multipart_files)} entries for {len(files)} files, skipping")
        return
    completed: list[CompletedMultipartSaveFile] = []
    failed: list[CompleteSaveFile] = []
    for file_info, multipart_info in zip(files, multipart_files):
        upload_id = multipart_info["uploadId"]
        parts = multipart_info["parts"]
        file_key = save_tracker_key(file_info)
        _track_file(tracker, file_key, file_info["source_path"], file_info["size"])
        completed_parts = _upload_multipart_parts(session, file_info, parts, part_size, on_part_progress)
        if completed_parts is None:
            console.warning(f"Multipart upload failed for {file_info['path']}, skipping")
            failed.append({"path": file_info["path"], "state": "FAILED"})
            continue
        if tracker is not None:
            tracker.finish_file(file_key)
        completed.append(
            {
                "path": file_info["path"],
                "parts": completed_parts,
                "uploadId": upload_id,
            }
        )
    if completed:
        with safe.block(message="Failed to complete multipart save upload"):
            complete_multipart_save(ctx.experiment_id, files=completed)
        console.info(f"Synced {len(completed)} save file(s)")
    if failed:
        with safe.block(message="Failed to report multipart save failures"):
            complete_save_files(ctx.experiment_id, files=failed)


def _upload_multipart_parts(
    session: SessionWithRetry,
    file_info: SaveFileEntry,
    parts: Sequence[MultipartPart],
    part_size: int,
    on_part_progress: Optional[PartProgressReporter],
) -> Optional[list[CompletedPart]]:
    """上传单个大文件的所有分片，返回 CompletedPart 列表。"""
    completed_parts: list[CompletedPart] = []
    source_path = file_info["source_path"]
    file_key = save_tracker_key(file_info)

    def upload_part(part_info: MultipartPart) -> CompletedPart:
        part_number = part_info["partNumber"]
        url = part_info["url"]
        offset = (part_number - 1) * part_size
        part_key = f"{file_key}:part:{part_number}"
        part_length = min(part_size, max(file_info["size"] - offset, 0))

        def on_read(current_bytes):
            if on_part_progress is not None:
                on_part_progress(file_key, part_key, current_bytes)

        try:
            with open(source_path, "rb") as f:
                wrapped_data = ProgressFileWrapper(f, on_read, part_length, offset=offset, size=part_length)
                resp = session.put(url, data=wrapped_data, headers={"Content-Type": "application/octet-stream"})
            resp.raise_for_status()
            if on_part_progress is not None:
                on_part_progress(file_key, part_key, part_length, final=True)
            etag = resp.headers.get("ETag", "").strip('"')
            return {"partNumber": part_number, "etag": etag}
        except Exception:
            if on_part_progress is not None:
                on_part_progress(file_key, part_key, 0)
            raise

    failed = False
    with SafeThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(upload_part, p) for p in parts]
        for future in futures:
            result: Optional[CompletedPart] = None
            with safe.block(message=f"Multipart part failed for {file_info['path']}"):
                result = future.result()
            if result is None:
                failed = True
                break
            completed_parts.append(result)
    if failed:
        console.warning(f"Multipart upload aborted for {file_info['path']}, reporting failure to server")
        return None
    return completed_parts


__all__ = [
    "PartProgressReporter",
    "put_media_resources",
    "upload_save_files",
]
