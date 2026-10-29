"""
@author: caddiesnew
@file: __init__.py
@time: 2026/3/19
@description: Record HTTP 上传接口
"""

from __future__ import annotations

import json
import math
import threading
from collections.abc import Callable, Iterable, Sequence
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Union, cast

import requests
import yaml

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.sdk.internal.core_python.api.upload import (
    upload_columns,
    upload_conda,
    upload_config,
    upload_histogram,
    upload_log,
    upload_media,
    upload_metadata,
    upload_requirements,
    upload_scalar,
)
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.pkg import adapter, client, console
from swanlab.sdk.typings.core_python.api.upload import (
    UploadColumn,
    UploadHistogram,
    UploadHistogramBatch,
    UploadLog,
    UploadLogBatch,
    UploadMediaBatch,
    UploadScalar,
    UploadScalarBatch,
)
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult

from .media import MediaPayload, MediaRecordTask, collect_media_tasks, put_media_resources
from .save import resolve_save_source, upload_save_files
from .utils import encode_column, is_retryable_error, transport_failure

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
    from swanlab.sdk.internal.pkg.client.session import SessionWithRetry


class HttpRecordSender:
    """Record HTTP 上传接口，按类型分发请求并返回 UploadResult。"""

    def __init__(self, ctx: CoreContext) -> None:
        self._ctx = ctx
        self._username = ctx.username
        self._project = ctx.project
        self._project_id = ctx.project_id
        self._experiment_id = ctx.experiment_id
        self._tracker: Optional[UploadTracker] = None
        self._completed_upload_keys: set[str] = set()
        self._completed_upload_keys_lock = threading.Lock()
        self._buffer_session: Optional[SessionWithRetry] = None
        self._upload_handlers: dict[str, Callable[[Sequence[Record]], UploadResult]] = {
            "column": self.upload_column,
            "scalar": self.upload_scalar,
            "media": self.upload_media,
            "histogram": self.upload_histogram,
            "log": self.upload_log,
            "save": self.upload_save,
        }

    def set_tracker(self, tracker: UploadTracker) -> None:
        self._tracker = tracker

    def close(self) -> None:
        """关闭对象存储上传 session。"""
        if self._buffer_session is not None:
            self._buffer_session.close()
            self._buffer_session = None

    def _ensure_session(self) -> SessionWithRetry:
        """获取对象存储上传 session（default_retry=0，由 transport 层控制重试）。"""
        if self._buffer_session is None:
            self._buffer_session = client.session.create(default_retry=0)
        return self._buffer_session

    def _mark_upload_progress(self, file_key: str, upload_key: str, current_bytes: int, *, final: bool = False) -> None:
        """上报分片字节进度，跳过完成的 upload_key。"""
        with self._completed_upload_keys_lock:
            if upload_key in self._completed_upload_keys:
                return
            if final:
                self._completed_upload_keys.add(upload_key)
        if self._tracker is not None:
            self._tracker.update_file_progress(file_key, upload_key, current_bytes)

    def _blocked_failure(self, records: Sequence[Record], message: str) -> UploadFailure:
        """构造本地数据校验失败结果。"""
        return UploadFailure(message=message, records=list(records))

    def upload(self, record_type: str, records: Sequence[Record]) -> UploadResult:
        """按 record_type 分发上传，返回 UploadResult。"""
        if not records:
            return UploadResult()
        handler = self._upload_handlers.get(record_type)
        if handler is None:
            console.warning(
                f"No upload handler registered for record type {record_type!r}; skipping {len(records)} record(s)."
            )
            return UploadResult(
                confirmed_records=list(records),
                failures=[self._blocked_failure(records, f"no upload handler for record type {record_type!r}")],
            )
        try:
            return handler(records)
        except requests.RequestException as e:
            retryable = is_retryable_error(e)
            if not retryable:
                console.warning(f"Failed to upload {record_type} records, skipping; error: {e}")
            return UploadResult(
                confirmed_records=[] if retryable else list(records),
                remaining_records=list(records) if retryable else [],
                failures=[transport_failure(list(records), e)],
            )

    def upload_column(self, records: Sequence[Record], batch_size: int = 3000) -> UploadResult:
        """分批上传列定义。"""
        # 先完成全部编码再发请求：编码失败必须发生在首个请求之前，避免重发已成功的批次
        try:
            encoded: List[UploadColumn] = [encode_column(record.column) for record in records]
        except KeyError as e:
            return UploadResult(
                confirmed_records=list(records),
                failures=[self._blocked_failure(records, f"failed to encode column: {e}")],
            )
        confirmed: List[Record] = []
        failures: List[UploadFailure] = []
        for start in range(0, len(encoded), batch_size):
            batch = list(records[start : start + batch_size])
            columns = encoded[start : start + batch_size]
            try:
                upload_columns(self._username, self._project, columns={"series": columns})
            except requests.RequestException as e:
                failure = transport_failure(batch, e)
                failures.append(failure)
                if is_retryable_error(e):
                    return UploadResult(
                        confirmed_records=confirmed,
                        remaining_records=list(records[start:]),
                        failures=failures,
                    )
                # 确定性错误：剩余批次同样会失败，一并终结不再重复提交
                confirmed.extend(records[start:])
                return UploadResult(confirmed_records=confirmed, failures=failures)
            confirmed.extend(batch)
        return UploadResult(confirmed_records=confirmed, failures=failures)

    def upload_scalar(self, records: Sequence[Record]) -> UploadResult:
        metrics: UploadScalarBatch = []
        confirmed: List[Record] = []
        blocked: List[Record] = []
        for record in records:
            if not record.HasField("scalar") or not record.scalar.HasField("timestamp"):
                blocked.append(record)
                continue
            scalar_record = record.scalar
            create_time = scalar_record.timestamp.ToJsonString()
            data = scalar_record.value.number if math.isfinite(scalar_record.value.number) else "nan"
            metric: UploadScalar = {
                "key": scalar_record.key,
                "index": scalar_record.step,
                "data": cast(Union[float, Literal["nan"]], data),
                "create_time": create_time,
            }
            metrics.append(metric)
            confirmed.append(record)
        console.debug(f"HTTP upload: upload_scalar with metrics count: {len(metrics)}")
        if not metrics:
            return UploadResult(
                confirmed_records=list(records),
                failures=[self._blocked_failure(blocked or list(records), "scalar record without value or timestamp")],
            )
        try:
            upload_scalar(self._project_id, self._experiment_id, metrics=metrics)
        except requests.RequestException as e:
            failures = [transport_failure(confirmed, e)]
            if blocked:
                failures.append(self._blocked_failure(blocked, "scalar record without value or timestamp"))
            retryable = is_retryable_error(e)
            return UploadResult(
                confirmed_records=blocked + ([] if retryable else confirmed),
                remaining_records=confirmed if retryable else [],
                failures=failures,
            )
        return UploadResult(
            confirmed_records=confirmed + blocked,
            failures=[self._blocked_failure(blocked, "scalar record without value or timestamp")] if blocked else [],
        )

    def upload_log(self, records: Sequence[Record]) -> UploadResult:
        console.debug(f"Preparing log HTTP request: records={len(records)}", write_to_tty=False)
        metrics: UploadLogBatch = []
        confirmed: List[Record] = []
        blocked: List[Record] = []
        for record in records:
            if not record.HasField("log") or not record.log.HasField("timestamp"):
                blocked.append(record)
                continue
            log_record = record.log
            create_time = log_record.timestamp.ToJsonString()
            metric: UploadLog = {
                "level": adapter.level[log_record.level],
                "epoch": log_record.epoch,
                "message": log_record.line,
                "create_time": create_time,
            }
            metrics.append(metric)
            confirmed.append(record)
        console.debug(f"Sending log HTTP request: metrics={len(metrics)}", write_to_tty=False)
        if not metrics:
            return UploadResult(
                confirmed_records=list(records),
                failures=[self._blocked_failure(blocked or list(records), "log record without line or timestamp")],
            )
        try:
            upload_log(self._project_id, self._experiment_id, metrics=metrics)
        except requests.RequestException as e:
            failures = [transport_failure(confirmed, e)]
            if blocked:
                failures.append(self._blocked_failure(blocked, "log record without line or timestamp"))
            retryable = is_retryable_error(e)
            return UploadResult(
                confirmed_records=blocked + ([] if retryable else confirmed),
                remaining_records=confirmed if retryable else [],
                failures=failures,
            )
        console.debug(f"Log HTTP request completed: metrics={len(metrics)}", write_to_tty=False)
        return UploadResult(
            confirmed_records=confirmed + blocked,
            failures=[self._blocked_failure(blocked, "log record without line or timestamp")] if blocked else [],
        )

    def upload_histogram(self, records: Sequence[Record]) -> UploadResult:
        """上传直方图指标，data 携带每桶计数与桶边界两组数组。"""
        metrics: UploadHistogramBatch = []
        confirmed: List[Record] = []
        blocked: List[Record] = []
        for record in records:
            if not record.HasField("histogram") or not record.histogram.HasField("timestamp"):
                blocked.append(record)
                continue
            histogram_record = record.histogram
            metric: UploadHistogram = {
                "key": histogram_record.key,
                # Histogram 未记录训练轮次，显式使用 0，避免 House 转发 null。
                "epoch": 0,
                "index": histogram_record.step,
                "data": {
                    "_type": "histogram",
                    "values": list(histogram_record.value.values),
                    "bins": list(histogram_record.value.bins),
                },
                "create_time": histogram_record.timestamp.ToJsonString(),
            }
            metrics.append(metric)
            confirmed.append(record)
        console.debug(f"HTTP upload: upload_histogram with metrics count: {len(metrics)}")
        if not metrics:
            return UploadResult(
                confirmed_records=list(records),
                failures=[
                    self._blocked_failure(blocked or list(records), "histogram record without value or timestamp")
                ],
            )
        try:
            upload_histogram(self._project_id, self._experiment_id, metrics=metrics)
        except requests.RequestException as e:
            failures = [transport_failure(confirmed, e)]
            if blocked:
                failures.append(self._blocked_failure(blocked, "histogram record without value or timestamp"))
            retryable = is_retryable_error(e)
            return UploadResult(
                confirmed_records=blocked + ([] if retryable else confirmed),
                remaining_records=confirmed if retryable else [],
                failures=failures,
            )
        return UploadResult(
            confirmed_records=confirmed + blocked,
            failures=[self._blocked_failure(blocked, "histogram record without value or timestamp")] if blocked else [],
        )

    def upload_media(self, records: Sequence[Record]) -> UploadResult:
        """收集媒体上传任务，上传至对象存储后提交媒体指标。"""
        tasks = collect_media_tasks(self._ctx, records)
        payloads: List[MediaPayload] = []
        blocked_records: List[Record] = []
        blocked_failures: List[UploadFailure] = []
        for task in tasks:
            payloads.extend(task.payloads)
            if task.failure is not None:
                blocked_records.append(task.record)
                blocked_failures.append(task.failure)

        if not payloads:
            return UploadResult(
                confirmed_records=list(records),
                failures=blocked_failures or [self._blocked_failure(records, "no media payload collected")],
            )
        uploadable_records = [task.record for task in tasks if task.metric is not None]
        try:
            put_failures = put_media_resources(
                self._ensure_session(), self._experiment_id, payloads, tracker=self._tracker
            )
        except requests.RequestException as e:
            # 取预签名 URL 失败属于整批失败：网络/5xx 保留重试，确定性错误整批终结
            retryable = is_retryable_error(e)
            if not retryable:
                self._fail_task_files(tasks)
            return UploadResult(
                confirmed_records=blocked_records + ([] if retryable else uploadable_records),
                remaining_records=uploadable_records if retryable else [],
                failures=blocked_failures + [transport_failure(uploadable_records, e)],
            )
        except OSError as e:
            self._fail_task_files(tasks)
            return UploadResult(
                confirmed_records=list(records),
                failures=blocked_failures
                + [self._blocked_failure([task.record for task in tasks], f"{type(e).__name__}: {e}")],
            )

        # 对象存储 PUT 错误逐文件存在：按任务聚合，结果不随文件顺序变化
        errors_by_task: Dict[int, List[BaseException]] = {}
        for payload_index, error in put_failures.items():
            errors_by_task.setdefault(payloads[payload_index].task_index, []).append(error)

        confirmed: List[Record] = list(blocked_records)
        remaining: List[Record] = []
        failures: List[UploadFailure] = list(blocked_failures)
        ready: List[MediaRecordTask] = []
        for task in tasks:
            if task.metric is None:
                continue
            errors = errors_by_task.get(task.index)
            if errors is None:
                ready.append(task)
                continue
            for error in errors:
                failures.append(transport_failure([task.record], error))
            if any(not isinstance(error, requests.RequestException) for error in errors):
                # 本地文件不可读或未知异常：重试不会恢复，终结该记录
                confirmed.append(task.record)
                self._fail_task_files([task])
            else:
                # 对象存储问题：保留记录重试，下一轮重新取预签名 URL
                remaining.append(task.record)

        if ready:
            ready_records = [task.record for task in ready]
            # ready 中的任务一定有 metric
            metrics = cast(UploadMediaBatch, [task.metric for task in ready])
            try:
                upload_media(self._project_id, self._experiment_id, metrics=metrics)
            except requests.RequestException as e:
                retryable = is_retryable_error(e)
                failures.append(transport_failure(ready_records, e))
                if retryable:
                    remaining.extend(ready_records)
                else:
                    confirmed.extend(ready_records)
                    self._fail_task_files(ready)
            else:
                self._finish_task_files(ready)
                confirmed.extend(ready_records)

        return UploadResult(confirmed_records=confirmed, remaining_records=remaining, failures=failures)

    def _fail_task_files(self, tasks: Iterable[MediaRecordTask]) -> None:
        """关闭已终结媒体任务的文件进度行。"""
        if self._tracker is None:
            return
        for task in tasks:
            for payload in task.payloads:
                self._tracker.fail_file(payload.tracker_key)

    def _finish_task_files(self, tasks: Iterable[MediaRecordTask]) -> None:
        """标记已成功上传媒体任务的文件进度行。"""
        if self._tracker is None:
            return
        for task in tasks:
            for payload in task.payloads:
                self._tracker.finish_file(payload.tracker_key)

    def upload_save(self, records: Sequence[Record]) -> UploadResult:
        """分发内部保存文件与用户保存文件（CUSTOM）上传。"""
        confirmed: List[Record] = []
        remaining: List[Record] = []
        failures: List[UploadFailure] = []
        custom_records: List[Record] = []
        for record in records:
            if not record.HasField("save"):
                failures.append(self._blocked_failure([record], "record without save payload"))
                confirmed.append(record)
                continue
            save = record.save
            if save.type == SaveType.SAVE_TYPE_CUSTOM:
                if save.HasField("payload"):
                    console.warning(f"CUSTOM save must not carry payload, skipping: {save.name or save.source_path}")
                    failures.append(self._blocked_failure([record], "CUSTOM save must not carry payload"))
                    confirmed.append(record)
                    continue
                custom_records.append(record)
                continue
            result = self.upload_internal_save(record)
            confirmed.extend(result.confirmed_records)
            remaining.extend(result.remaining_records)
            failures.extend(result.failures)
        if custom_records:
            result = upload_save_files(
                self._ctx,
                self._ensure_session(),
                custom_records,
                tracker=self._tracker,
                on_part_progress=self._mark_upload_progress,
            )
            confirmed.extend(result.confirmed_records)
            remaining.extend(result.remaining_records)
            failures.extend(result.failures)
        else:
            console.debug("No user save records to upload after filtering; all saves are internal files.")
        return UploadResult(confirmed_records=confirmed, remaining_records=remaining, failures=failures)

    def upload_internal_save(self, record: Record) -> UploadResult:
        """读取、解码并上传单条内部 save 记录。"""
        save = record.save
        data = self.read_internal_save(save)
        if data is None:
            return UploadResult(
                confirmed_records=[record],
                failures=[self._blocked_failure([record], f"internal save content unavailable: {save.source_path}")],
            )
        content = self.decode_internal_save(save.type, data)
        if content is None:
            return UploadResult(
                confirmed_records=[record],
                failures=[self._blocked_failure([record], f"internal save payload unparseable: {save.type}")],
            )
        try:
            if save.type == SaveType.SAVE_TYPE_METADATA:
                upload_metadata(self._username, self._project, self._experiment_id, content=content)
            elif save.type == SaveType.SAVE_TYPE_REQUIREMENTS:
                upload_requirements(self._username, self._project, self._experiment_id, content=content)
            elif save.type == SaveType.SAVE_TYPE_CONDA:
                upload_conda(self._username, self._project, self._experiment_id, content=content)
            elif save.type == SaveType.SAVE_TYPE_CONFIG:
                upload_config(self._username, self._project, self._experiment_id, content=content)
        except requests.RequestException as e:
            retryable = is_retryable_error(e)
            return UploadResult(
                confirmed_records=[] if retryable else [record],
                remaining_records=[record] if retryable else [],
                failures=[transport_failure([record], e)],
            )
        return UploadResult(confirmed_records=[record])

    def read_internal_save(self, save: SaveRecord) -> Optional[bytes]:
        """读取内部 save 字节内容（payload 或本地文件）。"""
        if save.HasField("payload"):
            return save.payload
        if self._ctx.config.skip_store:
            console.warning(f"Internal save payload missing with skip_store enabled, skipping: {save.type}")
            return None
        source = resolve_save_source(self._ctx, save)
        if source is None:
            console.warning(f"Save file not found, skipping: {save.source_path}")
            return None
        try:
            return source.read_bytes()
        except OSError as e:
            console.warning(f"Failed to read save file, skipping; file kept at {source}: {e}")
            return None

    def decode_internal_save(self, save_type: SaveType, data: bytes) -> Optional[Any]:
        """按 SaveType 解码内部 save 字节内容。"""
        try:
            if save_type == SaveType.SAVE_TYPE_METADATA:
                content = json.loads(data.decode("utf-8"))
                return content if isinstance(content, dict) else None
            if save_type in (SaveType.SAVE_TYPE_REQUIREMENTS, SaveType.SAVE_TYPE_CONDA):
                text = data.decode("utf-8")
                return text or None
            if save_type == SaveType.SAVE_TYPE_CONFIG:
                content = yaml.safe_load(data.decode("utf-8"))
                return content if isinstance(content, dict) else None
        except (UnicodeDecodeError, json.JSONDecodeError, yaml.YAMLError) as e:
            console.warning(f"Failed to parse internal save payload, skipping: type={save_type}, error={e}")
            return None
        console.warning(f"Unknown internal save type {save_type}, skipping")
        return None


__all__ = [
    "HttpRecordSender",
]
