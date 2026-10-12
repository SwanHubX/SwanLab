"""
@author: caddiesnew
@file: files.py
@time: 2026/4/16
@description: 媒体与 CUSTOM save 的对象存储传输：presigned PUT 并发、分片上传、字节进度上报

本模块是自由函数集合，session/tracker 等依赖全部显式传参，不持有模块级状态；
Record 编码与收集见 utils.py，record 路由与业务指标提交见 sender/__init__.py。
错误不在此吞掉：网络与业务异常转为结果中的 UploadFailure 并映射回原 Record；
文件完成确认（finish_file）一律在对应后端确认阶段（complete / 媒体指标提交）成功后执行。
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, List, Optional, Protocol, Tuple

import requests

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
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.client.session import SessionWithRetry
from swanlab.sdk.internal.pkg.executor import SafeThreadPoolExecutor
from swanlab.sdk.typings.core_python.api.save import (
    CompletedPart,
    CompleteSaveFile,
    MultipartPart,
    PrepareMultipartSaveFile,
    PrepareSaveFile,
    SaveFileEntry,
)
from swanlab.sdk.typings.core_python.api.upload import UploadResource
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult

from .utils import MediaPayload, compute_md5, resolve_save_source, save_tracker_key, transport_failure

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
    """媒体载荷批量传输到对象存储（presigned PUT 并发），并注册字节进度显示行。

    传输失败向上抛出，由调用方转为失败结果；文件完成确认由调用方在媒体指标
    提交成功后执行，本函数只负责传输字节进度。
    """
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


@dataclass
class _SaveFileTask:
    """单个 CUSTOM save 文件的上传任务：API DTO 与原 Record、tracker key 的关联。"""

    record: Record
    entry: SaveFileEntry
    tracker_key: str


def upload_save_files(
    ctx: CoreContext,
    session: SessionWithRetry,
    records: Sequence[Record],
    *,
    tracker: Optional[UploadTracker] = None,
    on_part_progress: Optional[PartProgressReporter] = None,
) -> UploadResult:
    """CUSTOM save 文件上传编排：解析来源 → 校验大小 → 按 save_batch 分批 → 小/大文件分流。

    调用方保证 records 均为 CUSTOM 类型且不携带 payload。
    返回记录级结果：complete 确认成功的文件进入 confirmed；来源缺失、超限、协议不匹配、
    传输或 complete 失败的 Record 进入 remaining 并携带失败原因，等待整流程重试
    （重跑会重新 prepare——后端幂等并重置 state、发新地址；同 key 重复 PUT 为覆盖写）。
    """
    config = ctx.config
    confirmed: List[Record] = []
    remaining: List[Record] = []
    failures: List[UploadFailure] = []

    # 1. 收集合法文件（不读内容，避免大批次内存爆炸）
    pending: List[Tuple[str, int, str, Record]] = []  # (source_path, size, name, record)
    for record in records:
        if not record.HasField("save"):
            failures.append(UploadFailure(source="local_data", message="record without save payload", records=[record]))
            remaining.append(record)
            continue
        save = record.save
        source_ref_path = resolve_save_source(ctx, save)
        if source_ref_path is None:
            console.warning(f"Save file not found: {save.source_path}")
            failures.append(
                UploadFailure(source="local_data", message=f"save file not found: {save.source_path}", records=[record])
            )
            remaining.append(record)
            continue
        try:
            size = source_ref_path.stat().st_size
        except OSError as e:
            console.warning(f"Failed to stat save file, retained for retry: {source_ref_path}: {e}")
            failures.append(
                UploadFailure(
                    source="local_data",
                    message=f"save file unreadable: {source_ref_path}",
                    records=[record],
                    exception=e,
                )
            )
            remaining.append(record)
            continue
        if size > config.save_size:
            console.warning(f"Save file exceeds size limit ({size} > {config.save_size}): {source_ref_path}")
            failures.append(
                UploadFailure(
                    source="local_data", message=f"save file exceeds size limit: {source_ref_path}", records=[record]
                )
            )
            remaining.append(record)
            continue
        # 用解析后的可读路径（可能是回退到 files 子目录的路径）作为实际读取路径
        # 保持 str() 原生路径形式，与下游 compute_md5 / upload_saves 的路径形式一致，
        # 避免 Windows 上 .as_posix() 与 str() 形式不一致导致映射查不到
        pending.append((str(source_ref_path), size, save.name, record))
    if not pending:
        return UploadResult(remaining_records=remaining, failures=failures)
    if config.save_batch <= 0:
        # 配置错误属确定性阻塞：整批保留并报告，不发起请求
        console.warning(f"Invalid save batch size ({config.save_batch}), retaining save upload")
        batch_records = [record for *_, record in pending]
        failures.append(
            UploadFailure(
                source="local_data", message=f"invalid save batch size {config.save_batch}", records=batch_records
            )
        )
        return UploadResult(remaining_records=remaining + batch_records, failures=failures)

    # 2. 按 save_batch 拆分文件列表，避免单次 prepare/complete 超过后端限制
    for index in range(0, len(pending), config.save_batch):
        result = _upload_save_batch(
            ctx,
            session,
            pending[index : index + config.save_batch],
            tracker=tracker,
            on_part_progress=on_part_progress,
        )
        confirmed.extend(result.confirmed_records)
        remaining.extend(result.remaining_records)
        failures.extend(result.failures)
    return UploadResult(confirmed_records=confirmed, remaining_records=remaining, failures=failures)


def _upload_save_batch(
    ctx: CoreContext,
    session: SessionWithRetry,
    pending: Sequence[Tuple[str, int, str, Record]],
    *,
    tracker: Optional[UploadTracker] = None,
    on_part_progress: Optional[PartProgressReporter] = None,
) -> UploadResult:
    """上传一个 save_batch 内的文件列表，返回记录级结果。"""
    config = ctx.config
    confirmed: List[Record] = []
    remaining: List[Record] = []
    failures: List[UploadFailure] = []
    tasks: List[_SaveFileTask] = []
    with SafeThreadPoolExecutor(max_workers=4) as executor:
        md5_futures = {source_path: executor.submit(compute_md5, source_path) for source_path, _, _, _ in pending}
        for source_path, size, name, record in pending:
            try:
                md5 = md5_futures[source_path].result()
            except Exception as e:
                console.warning(f"Save file MD5 computation failed, retained for retry: {source_path}: {e}")
                failures.append(
                    UploadFailure(
                        source="local_data",
                        message=f"save file md5 computation failed: {source_path}",
                        records=[record],
                        exception=e,
                    )
                )
                remaining.append(record)
                continue
            entry: SaveFileEntry = {
                "path": name,
                "source_path": source_path,
                "size": size,
                "md5": md5,
                "mime_type": guess_type(source_path),
            }
            tasks.append(_SaveFileTask(record=record, entry=entry, tracker_key=save_tracker_key(entry)))

    small_tasks = [task for task in tasks if task.entry["size"] < config.save_split]
    large_tasks = [task for task in tasks if task.entry["size"] >= config.save_split]
    if small_tasks:
        result = _upload_small_saves(ctx, session, small_tasks, tracker=tracker)
        confirmed.extend(result.confirmed_records)
        remaining.extend(result.remaining_records)
        failures.extend(result.failures)
    if large_tasks:
        result = _upload_large_saves(ctx, session, large_tasks, tracker=tracker, on_part_progress=on_part_progress)
        confirmed.extend(result.confirmed_records)
        remaining.extend(result.remaining_records)
        failures.extend(result.failures)
    return UploadResult(confirmed_records=confirmed, remaining_records=remaining, failures=failures)


def _upload_small_saves(
    ctx: CoreContext,
    session: SessionWithRetry,
    tasks: List[_SaveFileTask],
    *,
    tracker: Optional[UploadTracker] = None,
) -> UploadResult:
    """小文件：prepare → presigned PUT（并发）→ complete（仅确认上传成功的文件）。

    PUT 失败的文件不再向后端补报 FAILED（后端对 FAILED 会删除对象），保留在结果中
    等待整流程重试；complete 只包含 PUT 成功的文件，complete 成功后才确认并结束文件行。
    """
    prepare_items: List[PrepareSaveFile] = [
        {
            "path": task.entry["path"],
            "size": task.entry["size"],
            "md5": task.entry["md5"],
            "mimeType": task.entry["mime_type"],
        }
        for task in tasks
    ]
    try:
        resp = prepare_save_files(ctx.experiment_id, files=prepare_items)
    except requests.RequestException as e:
        return UploadResult(
            remaining_records=[task.record for task in tasks],
            failures=[
                transport_failure([task.record for task in tasks], e, source="swanlab_api", stage="request_address")
            ],
        )
    urls = resp["urls"]
    if len(urls) != len(tasks):
        return UploadResult(
            remaining_records=[task.record for task in tasks],
            failures=[
                UploadFailure(
                    source="swanlab_api",
                    stage="request_address",
                    message=f"save prepare returned {len(urls)} URLs for {len(tasks)} files",
                    records=[task.record for task in tasks],
                )
            ],
        )

    resources: List[UploadResource] = []
    for task, url in zip(tasks, urls):
        _track_file(tracker, task.tracker_key, task.entry["source_path"], task.entry["size"])
        resources.append(
            {
                "url": url,
                "source_path": task.entry["source_path"],
                "content_type": task.entry["mime_type"],
                "size": task.entry["size"],
                "tracker_key": task.tracker_key,
            }
        )
    put_failures = upload_saves(session, resources=resources, tracker=tracker)

    uploaded: List[_SaveFileTask] = []
    failures: List[UploadFailure] = []
    for index, task in enumerate(tasks):
        error = put_failures.get(index)
        if error is None:
            uploaded.append(task)
            continue
        failures.append(
            transport_failure(
                [task.record],
                error,
                source="object_storage",
                stage="transfer",
                message=f"save file upload failed: {task.entry['path']}",
            )
        )
    if not uploaded:
        return UploadResult(remaining_records=[task.record for task in tasks], failures=failures)

    complete_items: List[CompleteSaveFile] = [{"path": task.entry["path"], "state": "UPLOADED"} for task in uploaded]
    try:
        complete_save_files(ctx.experiment_id, files=complete_items)
    except requests.RequestException as e:
        # complete 失败：内容可能已在对象存储，但未获后端确认 → 全部保留等待整流程重试
        return UploadResult(
            remaining_records=[task.record for task in tasks],
            failures=failures
            + [transport_failure([task.record for task in uploaded], e, source="swanlab_api", stage="complete")],
        )
    for task in uploaded:
        if tracker is not None:
            tracker.finish_file(task.tracker_key)
    console.info(f"Synced {len(uploaded)} save file(s)")
    return UploadResult(
        confirmed_records=[task.record for task in uploaded],
        remaining_records=[task.record for index, task in enumerate(tasks) if index in put_failures],
        failures=failures,
    )


def _upload_large_saves(
    ctx: CoreContext,
    session: SessionWithRetry,
    tasks: List[_SaveFileTask],
    *,
    tracker: Optional[UploadTracker] = None,
    on_part_progress: Optional[PartProgressReporter] = None,
) -> UploadResult:
    """大文件：prepare-multipart → 分片 PUT（并发，收齐全部结果）→ 批量 complete-multipart。

    分片或 complete 失败的文件保留在结果中等待整流程重跑（重新 prepare 会生成新 uploadId，
    不复用旧会话）；分片全部成功的文件在 complete-multipart 成功后才确认并结束文件行。
    """
    part_size = ctx.config.save_part
    prepare_items: List[PrepareMultipartSaveFile] = []
    for task in tasks:
        part_count = math.ceil(task.entry["size"] / part_size)
        prepare_items.append(
            {
                "path": task.entry["path"],
                "size": task.entry["size"],
                "md5": task.entry["md5"],
                "mimeType": task.entry["mime_type"],
                "count": part_count,
            }
        )
    try:
        resp = prepare_multipart_save(ctx.experiment_id, files=prepare_items)
    except requests.RequestException as e:
        return UploadResult(
            remaining_records=[task.record for task in tasks],
            failures=[
                transport_failure([task.record for task in tasks], e, source="swanlab_api", stage="request_address")
            ],
        )
    multipart_files = resp["files"]
    if len(multipart_files) != len(tasks):
        return UploadResult(
            remaining_records=[task.record for task in tasks],
            failures=[
                UploadFailure(
                    source="swanlab_api",
                    stage="request_address",
                    message=f"multipart prepare returned {len(multipart_files)} entries for {len(tasks)} files",
                    records=[task.record for task in tasks],
                )
            ],
        )

    confirmed: List[Record] = []
    remaining: List[Record] = []
    failures: List[UploadFailure] = []
    ready: List[Tuple[_SaveFileTask, str, List[CompletedPart]]] = []
    for task, multipart_info in zip(tasks, multipart_files):
        upload_id = multipart_info["uploadId"]
        parts = multipart_info["parts"]
        _track_file(tracker, task.tracker_key, task.entry["source_path"], task.entry["size"])
        completed_parts, part_error = _upload_multipart_parts(session, task, parts, part_size, on_part_progress)
        if part_error is not None:
            failures.append(
                transport_failure(
                    [task.record],
                    part_error,
                    source="object_storage",
                    stage="transfer",
                    message=f"multipart upload failed: {task.entry['path']}",
                )
            )
            remaining.append(task.record)
            continue
        ready.append((task, upload_id, completed_parts))

    if not ready:
        return UploadResult(remaining_records=remaining, failures=failures)

    try:
        complete_multipart_save(
            ctx.experiment_id,
            files=[
                {"path": task.entry["path"], "parts": parts, "uploadId": upload_id} for task, upload_id, parts in ready
            ],
        )
    except requests.RequestException as e:
        # complete-multipart 失败：分片可能已在对象存储，但未获后端确认 → 保留等待整流程重跑
        return UploadResult(
            remaining_records=remaining + [task.record for task, _, _ in ready],
            failures=failures
            + [transport_failure([task.record for task, _, _ in ready], e, source="swanlab_api", stage="complete")],
        )
    for task, _, _ in ready:
        if tracker is not None:
            tracker.finish_file(task.tracker_key)
        confirmed.append(task.record)
    console.info(f"Synced {len(ready)} save file(s)")
    return UploadResult(confirmed_records=confirmed, remaining_records=remaining, failures=failures)


def _upload_multipart_parts(
    session: SessionWithRetry,
    task: _SaveFileTask,
    parts: Sequence[MultipartPart],
    part_size: int,
    on_part_progress: Optional[PartProgressReporter],
) -> Tuple[List[CompletedPart], Optional[BaseException]]:
    """上传单个大文件的所有分片，返回 (已完成分片, 首个错误)。

    收集全部并发结果后才返回，不因第一个异常丢掉其余分片的完成状态；
    错误非 None 时文件整体未完成，已完成分片仅用于字节进度展示。
    """
    completed_parts: List[CompletedPart] = []
    source_path = task.entry["source_path"]
    file_key = task.tracker_key
    first_error: Optional[BaseException] = None

    def upload_part(part_info: MultipartPart) -> CompletedPart:
        part_number = part_info["partNumber"]
        url = part_info["url"]
        offset = (part_number - 1) * part_size
        part_key = f"{file_key}:part:{part_number}"
        part_length = min(part_size, max(task.entry["size"] - offset, 0))

        def on_read(current_bytes: int) -> None:
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

    with SafeThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(upload_part, p) for p in parts]
        for future in futures:
            try:
                completed_parts.append(future.result())
            except Exception as e:
                if first_error is None:
                    first_error = e
    return completed_parts, first_error


__all__ = [
    "PartProgressReporter",
    "put_media_resources",
    "upload_save_files",
]
