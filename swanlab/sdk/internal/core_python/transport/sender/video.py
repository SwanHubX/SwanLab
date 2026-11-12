"""视频逐文件上传状态；由 Transport 的视频工作线程调用。"""

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Sequence, cast

import requests

from swanlab.sdk.internal.core_python.api.upload import upload_media
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.client.session import TimeoutHTTPAdapter
from swanlab.sdk.typings.core_python.api.upload import UploadMedia
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult

from .media import MediaPayload, MediaRecordTask, collect_media_tasks, put_media_resources
from .utils import is_retryable_error

if TYPE_CHECKING:
    from swanlab.proto.swanlab.record.v1.record_pb2 import Record
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker


@dataclass
class _VideoTask:
    task: MediaRecordTask
    status: list[str]
    retry_at: dict[int, float] = field(default_factory=dict)
    rounds: dict[int, int] = field(default_factory=dict)
    metric: Optional[UploadMedia] = None


class VideoUploader:
    def __init__(self, ctx: CoreContext, tracker: Optional[UploadTracker]):
        from swanlab.sdk.internal.core_python.transport.worker.utils import RetryPolicy

        self._ctx = ctx
        self._tracker = tracker
        self._retry = RetryPolicy()
        self._tasks: dict[int, _VideoTask] = {}
        self._uploaded: set[tuple[str, str, int]] = set()
        # PUT 不经过业务 API 会话，避免记录预签名 URL；超时与现有上传会话一致。
        self._session = requests.Session()
        for scheme in ("http://", "https://"):
            self._session.mount(scheme, TimeoutHTTPAdapter(timeout=60, max_retries=0))

    def close(self) -> None:
        self._session.close()
        self._tasks.clear()
        self._uploaded.clear()

    @staticmethod
    def _retryable(error: BaseException) -> bool:
        if not isinstance(error, requests.RequestException):
            return False
        if isinstance(
            error,
            (requests.exceptions.InvalidURL, requests.exceptions.InvalidSchema, requests.exceptions.MissingSchema),
        ):
            return False
        response = error.response
        if response is not None:
            if response.status_code in (408, 425, 429) or response.status_code >= 500:
                return True
            # 过期签名重新申请；永久 403 终止该文件。
            if response.status_code in (400, 401, 403):
                body = response.content[:8192].lower()
                return any(
                    code in body
                    for code in (b"expiredtoken", b"requestexpired", b"request has expired", b"signatureexpired")
                )
            if response.status_code >= 400:
                return False
        return is_retryable_error(error)

    def _failure(self, record: Record, error: BaseException) -> UploadFailure:
        response = getattr(error, "response", None)
        status = response.status_code if response is not None else None
        # 不保存原始异常，requests 异常可能携带签名 URL。
        message = f"Video upload failed: {type(error).__name__}"
        if status is not None:
            message += f" (HTTP {status})"
        return UploadFailure(message=message, records=[record], status=status)

    @staticmethod
    def _cache_key(payload: MediaPayload, task: MediaRecordTask) -> tuple[str, str, int]:
        # collect_media_tasks 可能已过滤缺失文件，按路径精确找到原 item。
        filename = Path(payload.remote_path).name
        sha256 = next(
            (item.sha256 for item in task.record.media.value.items if Path(item.filename).name == filename), ""
        )
        return payload.remote_path, sha256, payload.size

    @staticmethod
    def _validate(payload: MediaPayload, task: MediaRecordTask) -> None:
        from swanlab.sdk.internal.run.transforms.video.source import CHUNK_SIZE, MIME_TYPES, check_size, validate_header

        filename = Path(payload.remote_path).name
        if not filename or any(item.filename != Path(item.filename).name for item in task.record.media.value.items):
            raise ValueError("Invalid video filename")
        item = next(item for item in task.record.media.value.items if item.filename == filename)
        fmt = Path(filename).suffix.lower().lstrip(".")
        if fmt not in MIME_TYPES:
            raise ValueError("Unsupported video format")
        payload.content_type = MIME_TYPES[fmt]
        check_size(payload.size)
        stream = payload.buffer
        opened = isinstance(stream, (str, Path))
        file = open(stream, "rb") if opened else stream
        position = file.tell()
        try:
            file.seek(0)
            digest = hashlib.sha256()
            size = 0
            header = file.read(65536)
            validate_header(header, fmt)
            digest.update(header)
            size += len(header)
            for chunk in iter(lambda: file.read(CHUNK_SIZE), b""):
                size += len(chunk)
                check_size(size)
                digest.update(chunk)
            check_size(size)
            if size != item.size or digest.hexdigest() != item.sha256:
                raise ValueError("Video file size or SHA-256 changed before upload")
        finally:
            if opened:
                file.close()
            else:
                file.seek(position)

    def delay(self) -> float:
        """存在就绪文件时立即调度，否则等待最早重试时间。"""
        now = time.monotonic()
        deadlines = []
        for state in self._tasks.values():
            if state.metric is not None:
                deadlines.append(state.retry_at.get(-1, now))
            else:
                deadlines.extend(
                    state.retry_at.get(i, now) for i, status in enumerate(state.status) if status == "pending"
                )
        return max(0.0, min(deadlines) - now) if deadlines else 0.0

    def _defer(self, state: _VideoTask, index: int) -> None:
        rounds = state.rounds.get(index, 0) + 1
        state.rounds[index] = rounds
        state.retry_at[index] = time.monotonic() + self._retry.schedule(rounds)

    def _close_progress(self, state: _VideoTask, *, failed: bool = False) -> None:
        if self._tracker is None:
            return
        close = self._tracker.fail_file if failed else self._tracker.finish_file
        for status, payload in zip(state.status, state.task.payloads):
            if status == "success":
                close(payload.tracker_key)

    def upload(self, records: Sequence[Record]) -> UploadResult:
        result = UploadResult()
        for record in records:
            if id(record) not in self._tasks:
                task = collect_media_tasks(self._ctx, [record])[0]
                missing = len(record.media.value.items) - len(task.payloads)
                if task.failure is not None:
                    console.warning(f"Video record '{record.media.key}' skipped: {task.failure.message}")
                elif missing:
                    console.warning(f"Video key '{record.media.key}': skipped {missing} unavailable file(s)")
                self._tasks[id(record)] = _VideoTask(task, ["pending"] * len(task.payloads))
        uploaded_one = False
        for record in records:
            state = self._tasks[id(record)]
            task = state.task
            if task.failure is not None:
                result.failures.append(task.failure)
                result.confirmed_records.append(record)
                continue
            for index, payload in enumerate(task.payloads):
                if state.status[index] != "pending" or state.retry_at.get(index, 0) > time.monotonic():
                    continue
                cache_key = self._cache_key(payload, task)
                if cache_key in self._uploaded:
                    state.status[index] = "success"
                    continue
                if uploaded_one:
                    continue
                try:
                    self._validate(payload, task)
                    if not isinstance(payload.buffer, (str, Path)):
                        payload.buffer.seek(0)
                    uploaded_one = True
                    errors = put_media_resources(
                        self._session, self._ctx.experiment_id, [payload], tracker=self._tracker
                    )
                    if errors:
                        raise errors[0]
                except Exception as error:
                    failure = self._failure(record, error)
                    result.failures.append(failure)
                    if self._retryable(error):
                        self._defer(state, index)
                    else:
                        state.status[index] = "failed"
                        console.warning(f"{failure.message}; skipped one video in '{record.media.key}'")
                        if self._tracker is not None:
                            self._tracker.fail_file(payload.tracker_key)
                else:
                    state.status[index] = "success"
                    self._uploaded.add(cache_key)
            if "pending" in state.status:
                result.remaining_records.append(record)
                continue
            if state.metric is None:
                good = [i for i, status in enumerate(state.status) if status == "success"]
                if not good or task.metric is None:
                    result.confirmed_records.append(record)
                    continue
                state.metric = cast(UploadMedia, dict(task.metric))
                state.metric["data"] = [task.metric["data"][i] for i in good]
                state.metric["more"] = [task.metric["more"][i] for i in good]
            if state.retry_at.get(-1, 0) > time.monotonic():
                result.remaining_records.append(record)
                continue
            try:
                assert state.metric is not None
                upload_media(self._ctx.project_id, self._ctx.experiment_id, metrics=[state.metric])
            except requests.RequestException as error:
                failure = self._failure(record, error)
                result.failures.append(failure)
                if is_retryable_error(error):
                    self._defer(state, -1)
                    result.remaining_records.append(record)
                else:
                    console.warning(f"{failure.message}; local video files retained")
                    self._close_progress(state, failed=True)
                    result.confirmed_records.append(record)
            else:
                self._close_progress(state)
                result.confirmed_records.append(record)
        for record in result.confirmed_records:
            self._tasks.pop(id(record), None)
        return result
