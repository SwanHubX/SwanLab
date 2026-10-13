"""
@author: caddiesnew
@file: save.py
@time: 2026/4/16
@description: CUSTOM save 的对象存储传输与 save 来源解析
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path, PureWindowsPath
from typing import TYPE_CHECKING, List, Optional, Tuple, cast

import requests

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.sdk.internal.core_python.api.save import (
    complete_multipart_save,
    complete_save_files,
    prepare_multipart_save,
    prepare_save_files,
)
from swanlab.sdk.internal.core_python.api.upload import upload_saves
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg.mime import guess_type
from swanlab.sdk.internal.core_python.utils import ProgressFileWrapper
from swanlab.sdk.internal.pkg import console, safe
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

from .utils import PartProgressReporter, compute_md5, is_retryable_error, track_file, transport_failure

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


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


def save_tracker_key(file_info: SaveFileEntry) -> str:
    """构造 save 文件在 tracker 中的注册 key。"""
    return f"{file_info['path']}:{file_info['size']}:{file_info['md5']}"


@dataclass
class _SaveFileTask:
    """单个 CUSTOM save 文件的上传任务。"""

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
    """解析、校验并分批上传 CUSTOM save 文件，返回 UploadResult。"""
    config = ctx.config
    confirmed: List[Record] = []
    remaining: List[Record] = []
    failures: List[UploadFailure] = []

    pending: List[Tuple[str, int, str, Record]] = []  # (source_path, size, name, record)
    for record in records:
        if not record.HasField("save"):
            failures.append(UploadFailure(source="local_data", message="record without save payload", records=[record]))
            confirmed.append(record)
            continue
        save = record.save
        source_ref_path = resolve_save_source(ctx, save)
        if source_ref_path is None:
            console.warning(f"Save file not found: {save.source_path}")
            failures.append(
                UploadFailure(source="local_data", message=f"save file not found: {save.source_path}", records=[record])
            )
            confirmed.append(record)
            continue
        try:
            size = source_ref_path.stat().st_size
        except OSError as e:
            console.warning(f"Failed to stat save file, skipping: {source_ref_path}: {e}")
            failures.append(
                UploadFailure(
                    source="local_data",
                    message=f"save file unreadable: {source_ref_path}",
                    records=[record],
                    exception=e,
                )
            )
            confirmed.append(record)
            continue
        if size > config.save_size:
            console.warning(f"Save file exceeds size limit ({size} > {config.save_size}): {source_ref_path}")
            failures.append(
                UploadFailure(
                    source="local_data", message=f"save file exceeds size limit: {source_ref_path}", records=[record]
                )
            )
            confirmed.append(record)
            continue
        pending.append((str(source_ref_path), size, save.name, record))
    if not pending:
        return UploadResult(confirmed_records=confirmed, remaining_records=remaining, failures=failures)
    if config.save_batch <= 0:
        console.warning(f"Invalid save batch size ({config.save_batch}), skipping save upload")
        batch_records = [record for *_, record in pending]
        failures.append(
            UploadFailure(
                source="local_data", message=f"invalid save batch size {config.save_batch}", records=batch_records
            )
        )
        return UploadResult(confirmed_records=confirmed + batch_records, remaining_records=remaining, failures=failures)

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
    """上传单个批次内的 save 文件。"""
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
                console.warning(f"Save file MD5 computation failed, skipping: {source_path}: {e}")
                failures.append(
                    UploadFailure(
                        source="local_data",
                        message=f"save file md5 computation failed: {source_path}",
                        records=[record],
                        exception=e,
                    )
                )
                confirmed.append(record)
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
    """上传小文件：prepare → 并发 PUT → complete。"""
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
        retryable = is_retryable_error(e)
        return UploadResult(
            confirmed_records=[] if retryable else [task.record for task in tasks],
            remaining_records=[task.record for task in tasks] if retryable else [],
            failures=[
                transport_failure([task.record for task in tasks], e, source="swanlab_api", stage="request_address")
            ],
        )
    urls = resp["urls"]
    if len(urls) != len(tasks):
        # 数量不匹配是确定性协议错误，终结记录而非重试
        return UploadResult(
            confirmed_records=[task.record for task in tasks],
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
        track_file(tracker, task.tracker_key, task.entry["source_path"], task.entry["size"])
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
    skipped: List[_SaveFileTask] = []
    failures: List[UploadFailure] = []
    for index, task in enumerate(tasks):
        error = put_failures.get(index)
        if error is None:
            uploaded.append(task)
            continue
        console.warning(f"Save file upload failed, skipping: {task.entry['path']}")
        failures.append(
            transport_failure(
                [task.record],
                error,
                source="object_storage",
                stage="transfer",
                message=f"save file upload failed: {task.entry['path']}",
            )
        )
        skipped.append(task)

    # 统一通知后端：上传成功的标记 UPLOADED，失败的标记 FAILED 由后端清理对象存储
    complete_items: List[CompleteSaveFile] = cast(
        List[CompleteSaveFile],
        [{"path": task.entry["path"], "state": "UPLOADED"} for task in uploaded]
        + [{"path": task.entry["path"], "state": "FAILED"} for task in skipped],
    )
    if complete_items:
        try:
            complete_save_files(ctx.experiment_id, files=complete_items)
        except requests.RequestException as e:
            retryable = is_retryable_error(e)
            if not retryable and tracker is not None:
                for task in uploaded + skipped:
                    tracker.fail_file(task.tracker_key)
            return UploadResult(
                confirmed_records=[] if retryable else [task.record for task in uploaded + skipped],
                remaining_records=[task.record for task in uploaded + skipped] if retryable else [],
                failures=failures
                + [
                    transport_failure(
                        [task.record for task in uploaded + skipped],
                        e,
                        source="swanlab_api",
                        stage="complete",
                    )
                ],
            )
    for task in uploaded:
        if tracker is not None:
            tracker.finish_file(task.tracker_key)
    for task in skipped:
        if tracker is not None:
            tracker.fail_file(task.tracker_key)
    if uploaded:
        console.info(f"Synced {len(uploaded)} save file(s)")
    return UploadResult(
        confirmed_records=[task.record for task in uploaded + skipped],
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
    """上传大文件：prepare-multipart → 分片 PUT → complete-multipart。"""
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
        retryable = is_retryable_error(e)
        return UploadResult(
            confirmed_records=[] if retryable else [task.record for task in tasks],
            remaining_records=[task.record for task in tasks] if retryable else [],
            failures=[
                transport_failure([task.record for task in tasks], e, source="swanlab_api", stage="request_address")
            ],
        )
    multipart_files = resp["files"]
    if len(multipart_files) != len(tasks):
        # 数量不匹配是确定性协议错误，终结记录而非重试
        return UploadResult(
            confirmed_records=[task.record for task in tasks],
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
    skipped: List[_SaveFileTask] = []
    failures: List[UploadFailure] = []
    ready: List[Tuple[_SaveFileTask, str, List[CompletedPart]]] = []
    for task, multipart_info in zip(tasks, multipart_files):
        upload_id = multipart_info["uploadId"]
        parts = multipart_info["parts"]
        track_file(tracker, task.tracker_key, task.entry["source_path"], task.entry["size"])
        completed_parts, part_error = _upload_multipart_parts(session, task, parts, part_size, on_part_progress)
        if part_error is not None:
            console.warning(f"Multipart upload failed, skipping: {task.entry['path']}")
            failures.append(
                transport_failure(
                    [task.record],
                    part_error,
                    source="object_storage",
                    stage="transfer",
                    message=f"multipart upload failed: {task.entry['path']}",
                )
            )
            skipped.append(task)
            continue
        ready.append((task, upload_id, completed_parts))

    # 向后端报告分片 PUT 失败的文件，触发后端清理对象存储
    if skipped:
        failed_items: List[CompleteSaveFile] = [{"path": task.entry["path"], "state": "FAILED"} for task in skipped]
        with safe.block(requests.RequestException, message=None):
            # 报告失败不阻塞流程，后端依赖超时清理兜底
            complete_save_files(ctx.experiment_id, files=failed_items)

    if not ready:
        if tracker is not None:
            for task in skipped:
                tracker.fail_file(task.tracker_key)
        return UploadResult(
            confirmed_records=[task.record for task in skipped],
            failures=failures,
        )

    try:
        complete_multipart_save(
            ctx.experiment_id,
            files=[
                {"path": task.entry["path"], "parts": parts, "uploadId": upload_id} for task, upload_id, parts in ready
            ],
        )
    except requests.RequestException as e:
        # complete-multipart 遇到网络或 5xx 时保留重试，确定性 4xx 终结。
        retryable = is_retryable_error(e)
        if not retryable and tracker is not None:
            for task in skipped + [task for task, _, _ in ready]:
                tracker.fail_file(task.tracker_key)
        return UploadResult(
            confirmed_records=[task.record for task in skipped]
            + ([] if retryable else [task.record for task, _, _ in ready]),
            remaining_records=[task.record for task, _, _ in ready] if retryable else [],
            failures=failures
            + [transport_failure([task.record for task, _, _ in ready], e, source="swanlab_api", stage="complete")],
        )
    for task in skipped:
        if tracker is not None:
            tracker.fail_file(task.tracker_key)
    for task, _, _ in ready:
        if tracker is not None:
            tracker.finish_file(task.tracker_key)
        confirmed.append(task.record)
    if confirmed:
        console.info(f"Synced {len(confirmed)} save file(s)")
    return UploadResult(
        confirmed_records=confirmed + [task.record for task in skipped],
        failures=failures,
    )


def _upload_multipart_parts(
    session: SessionWithRetry,
    task: _SaveFileTask,
    parts: Sequence[MultipartPart],
    part_size: int,
    on_part_progress: Optional[PartProgressReporter],
) -> Tuple[List[CompletedPart], Optional[BaseException]]:
    """并发上传大文件的分片，返回 (完成分片, 首个错误)。"""
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
    "resolve_save_source",
    "save_tracker_key",
    "upload_save_files",
]
