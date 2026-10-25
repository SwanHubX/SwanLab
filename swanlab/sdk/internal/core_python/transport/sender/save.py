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
    # 1. 优先使用记录中的原始绝对路径
    primary = Path(save.source_path)
    if primary.is_file():
        return primary
    # 2. skip_store 模式没有本 run 的 files 镜像，禁止回退到可能属于旧 run 的同名文件
    if ctx.config.skip_store:
        return None
    # 3. 按 save 类型在 run 目录的 files 子目录中回退定位
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


def _prepare_failure(tasks: Sequence[_SaveFileTask], error: requests.RequestException) -> UploadResult:
    """prepare 请求失败：可重试错误保留记录，其余终结。"""
    records = [task.record for task in tasks]
    retryable = is_retryable_error(error)
    return UploadResult(
        confirmed_records=[] if retryable else records,
        remaining_records=records if retryable else [],
        failures=[transport_failure(records, error)],
    )


def _prepare_count_mismatch(tasks: Sequence[_SaveFileTask], *, kind: str, received: int) -> UploadResult:
    """prepare 返回数量不匹配是协议错误，终结整批记录而非重试。"""
    records = [task.record for task in tasks]
    unit = "URLs" if kind == "save" else "entries"
    console.warning(f"{kind.title()} prepare returned {received} {unit} for {len(records)} files, skipping")
    return UploadResult(
        confirmed_records=records,
        failures=[
            UploadFailure(
                message=f"{kind} prepare returned {received} entries for {len(records)} files",
                records=records,
            )
        ],
    )


def _prepare_protocol_failure(tasks: Sequence[_SaveFileTask], *, kind: str, message: str) -> UploadResult:
    """prepare 响应字段缺失或结构错误时终结整批记录，避免进入无界重试。"""
    records = [task.record for task in tasks]
    console.warning(f"{kind.title()} prepare response invalid, skipping: {message}")
    return UploadResult(
        confirmed_records=records,
        failures=[UploadFailure(message=f"{kind} prepare response invalid: {message}", records=records)],
    )


def _response_list(response: object, key: str) -> Optional[list]:
    """读取响应中的列表字段，缺失或类型错误时返回 None。"""
    if not isinstance(response, dict):
        return None
    value = response.get(key)
    return value if isinstance(value, list) else None


def _validated_parts(value: object, expected_count: int) -> Tuple[Optional[List[MultipartPart]], str]:
    """校验分片列表的结构、数量与编号，返回 (分片列表, 错误描述)。"""
    # 1. 校验列表类型与分片数量
    if not isinstance(value, list):
        return None, "missing parts list"
    if len(value) != expected_count:
        return None, f"returned {len(value)} parts for {expected_count} expected"
    # 2. 逐项校验 partNumber/url 并重建类型安全的列表
    parts: List[MultipartPart] = []
    for entry in value:
        if not isinstance(entry, dict):
            return None, "part entry is not an object"
        part_number = entry.get("partNumber")
        url = entry.get("url")
        if not isinstance(part_number, int) or not isinstance(url, str) or not url:
            return None, "part entry missing partNumber or url"
        parts.append({"partNumber": part_number, "url": url})
    # 3. 校验编号连续为 1..N，避免缺片、重复或越界
    if sorted(part["partNumber"] for part in parts) != list(range(1, expected_count + 1)):
        return None, "part numbers are not consecutive from 1"
    return parts, ""


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

    # 1. 逐条解析源文件：缺失、不可 stat、超限的记录直接终结，其余进入待上传列表
    pending: List[Tuple[str, int, str, Record]] = []  # (source_path, size, name, record)
    for record in records:
        if not record.HasField("save"):
            failures.append(UploadFailure(message="record without save payload", records=[record]))
            confirmed.append(record)
            continue
        save = record.save
        source_ref_path = resolve_save_source(ctx, save)
        if source_ref_path is None:
            console.warning(f"Save file not found, skipping: {save.source_path}")
            failures.append(UploadFailure(message=f"save file not found: {save.source_path}", records=[record]))
            confirmed.append(record)
            continue
        try:
            size = source_ref_path.stat().st_size
        except OSError as e:
            console.warning(f"Failed to stat save file, skipping: {source_ref_path}: {e}")
            failures.append(
                UploadFailure(
                    message=f"save file unreadable: {source_ref_path}",
                    records=[record],
                    exception=e,
                )
            )
            confirmed.append(record)
            continue
        if size > config.save_size:
            console.warning(f"Save file exceeds size limit ({size} > {config.save_size}), skipping: {source_ref_path}")
            failures.append(UploadFailure(message=f"save file exceeds size limit: {source_ref_path}", records=[record]))
            confirmed.append(record)
            continue
        pending.append((str(source_ref_path), size, save.name, record))
    # 2. 没有可上传文件，或 save_batch 配置非法时，直接汇总返回
    if not pending:
        return UploadResult(confirmed_records=confirmed, remaining_records=remaining, failures=failures)
    if config.save_batch <= 0:
        console.warning(f"Invalid save batch size ({config.save_batch}), skipping save upload")
        batch_records = [record for *_, record in pending]
        failures.append(UploadFailure(message=f"invalid save batch size {config.save_batch}", records=batch_records))
        return UploadResult(confirmed_records=confirmed + batch_records, remaining_records=remaining, failures=failures)

    # 3. 按 save_batch 分批上传，逐批合并 confirmed/remaining/failures
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
    # 1. 并发计算 MD5，计算失败的记录直接终结
    with SafeThreadPoolExecutor(max_workers=4) as executor:
        md5_futures = {source_path: executor.submit(compute_md5, source_path) for source_path, _, _, _ in pending}
        for source_path, size, name, record in pending:
            try:
                md5 = md5_futures[source_path].result()
            except Exception as e:
                console.warning(f"Save file MD5 computation failed, skipping: {source_path}: {e}")
                failures.append(
                    UploadFailure(
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

    # 2. 按 save_split 拆成“小文件”和“分片上传”两组
    small_tasks = [task for task in tasks if task.entry["size"] < config.save_split]
    large_tasks = [task for task in tasks if task.entry["size"] >= config.save_split]
    # 3. 分别上传两组，并合并确认、待重试与失败结果
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
    # 1. prepare 申请预签名 URL，并校验响应字段、数量与 URL 合法性
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
        return _prepare_failure(tasks, e)
    urls = _response_list(resp, "urls")
    if urls is None:
        return _prepare_protocol_failure(tasks, kind="save", message="missing 'urls' list")
    if len(urls) != len(tasks):
        return _prepare_count_mismatch(tasks, kind="save", received=len(urls))
    if not all(isinstance(url, str) and url for url in urls):
        return _prepare_protocol_failure(tasks, kind="save", message="urls must be non-empty strings")
    urls = cast(List[str], urls)

    # 2. 并发 PUT 到对象存储，按错误类型分为成功、待重试与本地终结
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
    storage_failed: List[_SaveFileTask] = []
    local_failed: List[_SaveFileTask] = []
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
                message=f"save file upload failed: {task.entry['path']}",
            )
        )
        if isinstance(error, requests.RequestException):
            # 对象存储问题：保留记录重试，下一轮重新 prepare 取新链接
            storage_failed.append(task)
        else:
            # 本地读取失败或未知异常：终结该文件
            console.warning(f"Failed to upload save file, skipping: {task.entry['path']}: {error}")
            local_failed.append(task)

    # 3. 通知后端：上传成功的标记 UPLOADED，本地失败的标记 FAILED，待重试的文件保持 PENDING
    complete_items: List[CompleteSaveFile] = cast(
        List[CompleteSaveFile],
        [{"path": task.entry["path"], "state": "UPLOADED"} for task in uploaded]
        + [{"path": task.entry["path"], "state": "FAILED"} for task in local_failed],
    )
    if complete_items:
        try:
            complete_save_files(ctx.experiment_id, files=complete_items)
        except requests.RequestException as e:
            retryable = is_retryable_error(e)
            pending = uploaded + storage_failed
            if tracker is not None:
                # local_failed 已终结，无论 complete 是否保留重试都关闭其进度行
                for task in local_failed:
                    tracker.fail_file(task.tracker_key)
                if not retryable:
                    for task in pending:
                        tracker.fail_file(task.tracker_key)
            return UploadResult(
                confirmed_records=[task.record for task in local_failed]
                + ([] if retryable else [task.record for task in pending]),
                remaining_records=[task.record for task in pending] if retryable else [],
                failures=failures + [transport_failure([task.record for task in pending], e)],
            )
    for task in uploaded:
        if tracker is not None:
            tracker.finish_file(task.tracker_key)
    for task in local_failed:
        if tracker is not None:
            tracker.fail_file(task.tracker_key)
    if uploaded:
        console.info(f"Synced {len(uploaded)} save file(s)")
    return UploadResult(
        confirmed_records=[task.record for task in uploaded + local_failed],
        remaining_records=[task.record for task in storage_failed],
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
    # 1. prepare-multipart：计算分片数、申请 uploadId 与分片 URL，并校验响应结构与数量
    part_size = ctx.config.save_part
    expected_part_counts = [math.ceil(task.entry["size"] / part_size) for task in tasks]
    prepare_items: List[PrepareMultipartSaveFile] = []
    for task, part_count in zip(tasks, expected_part_counts):
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
        return _prepare_failure(tasks, e)
    multipart_files = _response_list(resp, "files")
    if multipart_files is None:
        return _prepare_protocol_failure(tasks, kind="multipart", message="missing 'files' list")
    if len(multipart_files) != len(tasks):
        return _prepare_count_mismatch(tasks, kind="multipart", received=len(multipart_files))

    validated: List[Tuple[_SaveFileTask, str, List[MultipartPart]]] = []
    for task, expected_count, multipart_info in zip(tasks, expected_part_counts, multipart_files):
        if not isinstance(multipart_info, dict):
            return _prepare_protocol_failure(tasks, kind="multipart", message="file entry is not an object")
        upload_id = multipart_info.get("uploadId")
        if not isinstance(upload_id, str) or not upload_id:
            return _prepare_protocol_failure(
                tasks, kind="multipart", message=f"invalid uploadId for {task.entry['path']!r}"
            )
        parts, error = _validated_parts(multipart_info.get("parts"), expected_count)
        if parts is None:
            return _prepare_protocol_failure(
                tasks, kind="multipart", message=f"invalid parts for {task.entry['path']!r}: {error}"
            )
        validated.append((task, upload_id, parts))

    confirmed: List[Record] = []
    local_failed: List[_SaveFileTask] = []
    storage_failed: List[_SaveFileTask] = []
    failures: List[UploadFailure] = []
    ready: List[Tuple[_SaveFileTask, str, List[CompletedPart]]] = []
    # 2. 并发上传分片，按错误类型分为 ready、待重试与本地终结
    for task, upload_id, parts in validated:
        track_file(tracker, task.tracker_key, task.entry["source_path"], task.entry["size"])
        completed_parts, part_error = _upload_multipart_parts(session, task, parts, part_size, on_part_progress)
        if part_error is not None:
            failures.append(
                transport_failure(
                    [task.record],
                    part_error,
                    message=f"multipart upload failed: {task.entry['path']}",
                )
            )
            if isinstance(part_error, requests.RequestException):
                # 对象存储问题：保留记录重试，下一轮重新 prepare 取新 uploadId
                storage_failed.append(task)
            else:
                # 本地读取失败或未知异常：终结该文件
                console.warning(f"Multipart upload failed for {task.entry['path']}, skipping: {part_error}")
                local_failed.append(task)
            continue
        ready.append((task, upload_id, completed_parts))

    # 3. 上报本地失败文件、complete-multipart，并关闭对应的进度行
    if local_failed:
        failed_items: List[CompleteSaveFile] = [
            {"path": task.entry["path"], "state": "FAILED"} for task in local_failed
        ]
        with safe.block(requests.RequestException, message=None):
            complete_save_files(ctx.experiment_id, files=failed_items)
        if tracker is not None:
            for task in local_failed:
                tracker.fail_file(task.tracker_key)

    if not ready:
        return UploadResult(
            confirmed_records=[task.record for task in local_failed],
            remaining_records=[task.record for task in storage_failed],
            failures=failures,
        )

    ready_tasks = [task for task, _, _ in ready]
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
        pending = ready_tasks + storage_failed
        if not retryable and tracker is not None:
            for task in pending:
                tracker.fail_file(task.tracker_key)
        return UploadResult(
            confirmed_records=[task.record for task in local_failed]
            + ([] if retryable else [task.record for task in pending]),
            remaining_records=[task.record for task in pending] if retryable else [],
            failures=failures + [transport_failure([task.record for task in pending], e)],
        )
    for task, _, _ in ready:
        if tracker is not None:
            tracker.finish_file(task.tracker_key)
        confirmed.append(task.record)
    if confirmed:
        console.info(f"Synced {len(confirmed)} save file(s)")
    return UploadResult(
        confirmed_records=confirmed + [task.record for task in local_failed],
        remaining_records=[task.record for task in storage_failed],
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

    # 1. 定义单个分片的上传任务：计算偏移、流式上传并返回 ETag
    def upload_part(part_info: MultipartPart) -> CompletedPart:
        part_number = part_info["partNumber"]
        url = part_info["url"]
        # 根据 partNumber 计算分片在文件中的偏移与长度
        offset = (part_number - 1) * part_size
        part_key = f"{file_key}:part:{part_number}"
        part_length = min(part_size, max(task.entry["size"] - offset, 0))

        def on_read(current_bytes: int) -> None:
            if on_part_progress is not None:
                on_part_progress(file_key, part_key, current_bytes)

        try:
            # 流式读取该区间并 PUT，同时上报字节进度
            with open(source_path, "rb") as f:
                wrapped_data = ProgressFileWrapper(f, on_read, part_length, offset=offset, size=part_length)
                resp = session.put(url, data=wrapped_data, headers={"Content-Type": "application/octet-stream"})
            resp.raise_for_status()
            if on_part_progress is not None:
                on_part_progress(file_key, part_key, part_length, final=True)
            etag = resp.headers.get("ETag", "").strip('"')
            return {"partNumber": part_number, "etag": etag}
        except Exception:
            # 失败时把当前分片进度归零，交给上层决定重试或终结
            if on_part_progress is not None:
                on_part_progress(file_key, part_key, 0)
            raise

    # 2. 并发提交所有分片，不因单个分片失败而提前中断其它分片
    with SafeThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(upload_part, p) for p in parts]
        for future in futures:
            try:
                completed_parts.append(future.result())
            except Exception as e:
                if first_error is None:
                    first_error = e
    # 3. 返回已完成分片与首个错误，调用方据此决定 complete 或重试
    return completed_parts, first_error


__all__ = [
    "resolve_save_source",
    "save_tracker_key",
    "upload_save_files",
]
