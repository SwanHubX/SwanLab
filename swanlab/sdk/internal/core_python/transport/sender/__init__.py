"""
@author: caddiesnew
@file: __init__.py
@time: 2026/3/19
@description: Record 与后端对接的唯一上传接口（HTTP 传输层）
"""

from __future__ import annotations

import json
import math
import threading
from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Any, List, Literal, Optional, Union, cast

import requests
import yaml

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.sdk.internal.core_python.api.upload import (
    upload_columns,
    upload_conda,
    upload_config,
    upload_log,
    upload_media,
    upload_metadata,
    upload_requirements,
    upload_scalar,
)
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.pkg import adapter, client, console
from swanlab.sdk.typings.core_python.api.upload import (
    UploadLog,
    UploadLogBatch,
    UploadScalar,
    UploadScalarBatch,
)
from swanlab.sdk.typings.core_python.transport import UploadFailure, UploadResult

from .files import put_media_resources, upload_save_files
from .utils import collect_media_payloads, encode_column, resolve_save_source, transport_failure

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
    from swanlab.sdk.internal.pkg.client.session import SessionWithRetry


class HttpRecordSender:
    """Record 与后端对接的真实上传接口。

    全部 upload_* 方法显式定义于本类，每个方法是一个完整短流程：
    编码/收集（utils.py）→ 对象存储字节传输（files.py）→ 业务指标提交（api 层）。
    所有 upload_* 方法返回 UploadResult：网络与业务错误转为结果内的 UploadFailure，
    未确认阶段的 Record 保留在 remaining，不因失败丢弃；未预期异常（编程错误）仍然上抛，
    由 worker 兜底。记录数统计由 Dispatch 依据 confirmed_records 推进，本类不推进。
    重试调度、错误分类不属于本类，由 worker 决定何时再次调用 upload。
    """

    def __init__(self, ctx: CoreContext) -> None:
        self._ctx = ctx
        self._username = ctx.username
        self._project = ctx.project
        self._project_id = ctx.project_id
        self._experiment_id = ctx.experiment_id
        self._tracker: Optional[UploadTracker] = None
        # 分片字节进度的已完成 upload-key 去重集合，跨重试生效（见 _mark_upload_progress）
        self._completed_upload_keys: set[str] = set()
        self._completed_upload_keys_lock = threading.Lock()
        # 资源上传 session，作用在Sender对象中复用TCP链接
        self._buffer_session: Optional[SessionWithRetry] = None
        self._upload_handlers: dict[str, Callable[[Sequence[Record]], UploadResult]] = {
            "column": self.upload_column,
            "scalar": self.upload_scalar,
            "media": self.upload_media,
            "log": self.upload_log,
            "save": self.upload_save,
        }

    def set_tracker(self, tracker: UploadTracker) -> None:
        self._tracker = tracker

    def _ensure_session(self) -> SessionWithRetry:
        """懒初始化资源上传 session（复用 TCP 连接）。

        default_retry=0：对象存储 PUT 与分片子线程单次执行（plan §6.3），
        重试由 transport 层的退避调度统一控制。
        """
        if self._buffer_session is None:
            self._buffer_session = client.session.create(default_retry=0)
        return self._buffer_session

    def _mark_upload_progress(self, file_key: str, upload_key: str, current_bytes: int, *, final: bool = False) -> None:
        """推送分片字节进度到 tracker，跳过已完成的上传事件 key（跨重试去重）。

        final=True 时将 upload_key 标记为已完成，后续对该 key 的调用（包括中间进度）将被跳过。
        """
        with self._completed_upload_keys_lock:
            if upload_key in self._completed_upload_keys:
                return
            if final:
                self._completed_upload_keys.add(upload_key)
        if self._tracker is not None:
            self._tracker.update_file_progress(file_key, upload_key, current_bytes)

    def _blocked_failure(self, records: Sequence[Record], message: str) -> UploadFailure:
        """本地数据不可执行的确定性失败：记录保留并报告，等待重试策略暂停调度。"""
        return UploadFailure(source="local_data", message=message, records=list(records))

    # ── 通用上传入口 ──

    def upload(self, record_type: str, records: Sequence[Record]) -> UploadResult:
        """通用上传入口。空输入返回空结果；未知类型整批保留并报告。"""
        if not records:
            return UploadResult()
        handler = self._upload_handlers.get(record_type)
        if handler is None:
            console.warning(f"No upload handler registered for record type {record_type!r}.")
            return UploadResult(
                remaining_records=list(records),
                failures=[self._blocked_failure(records, f"no upload handler for record type {record_type!r}")],
            )
        try:
            return handler(records)
        except requests.RequestException as e:
            # 兜底：handler 内部未分类的上传异常（含 ApiError），整批保留
            return UploadResult(
                remaining_records=list(records),
                failures=[transport_failure(list(records), e, source="swanlab_api")],
            )

    # ── 普通记录 ──

    def upload_column(self, records: Sequence[Record], batch_size: int = 3000) -> UploadResult:
        """列上传：分批提交，前序批次成功即确认（内层部分成功保留）。"""
        confirmed: List[Record] = []
        for start in range(0, len(records), batch_size):
            batch = list(records[start : start + batch_size])
            columns = [encode_column(record.column) for record in batch]
            try:
                upload_columns(self._username, self._project, columns={"series": columns})
            except requests.RequestException as e:
                failure = transport_failure(batch, e, source="swanlab_api", stage="submit_metrics")
                return UploadResult(
                    confirmed_records=confirmed,
                    remaining_records=list(records[start:]),
                    failures=[failure],
                )
            confirmed.extend(batch)
        return UploadResult(confirmed_records=confirmed)

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
                remaining_records=list(records),
                failures=[self._blocked_failure(blocked or list(records), "scalar record without value or timestamp")],
            )
        try:
            upload_scalar(self._project_id, self._experiment_id, metrics=metrics)
        except requests.RequestException as e:
            return UploadResult(
                remaining_records=list(records),
                failures=[transport_failure(list(records), e, source="swanlab_api", stage="submit_metrics")],
            )
        return UploadResult(
            confirmed_records=confirmed,
            remaining_records=blocked,
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
                remaining_records=list(records),
                failures=[self._blocked_failure(blocked or list(records), "log record without line or timestamp")],
            )
        try:
            upload_log(self._project_id, self._experiment_id, metrics=metrics)
        except requests.RequestException as e:
            return UploadResult(
                remaining_records=list(records),
                failures=[transport_failure(list(records), e, source="swanlab_api", stage="submit_metrics")],
            )
        console.debug(f"Log HTTP request completed: metrics={len(metrics)}", write_to_tty=False)
        return UploadResult(
            confirmed_records=confirmed,
            remaining_records=blocked,
            failures=[self._blocked_failure(blocked, "log record without line or timestamp")] if blocked else [],
        )

    # ── 媒体 ──

    def upload_media(self, records: Sequence[Record]) -> UploadResult:
        """媒体：收集载荷 → 对象存储 PUT → 提交媒体指标，两阶段均成功才算完成。

        文件行完成确认在指标提交成功后执行（plan §4.2）；未产出指标的 Record 保留并报告。
        """
        plan = collect_media_payloads(self._ctx, records)
        if not plan.payloads:
            # 全部 Record 均无可用载荷：保留并报告，不发起请求
            return UploadResult(
                remaining_records=list(records),
                failures=plan.blocked_failures or [self._blocked_failure(records, "no media payload collected")],
            )
        try:
            put_media_resources(self._ensure_session(), self._experiment_id, plan.payloads, tracker=self._tracker)
        except requests.RequestException as e:
            return UploadResult(
                remaining_records=list(records),
                failures=plan.blocked_failures
                + [transport_failure(list(records), e, source="object_storage", stage="transfer")],
            )
        try:
            upload_media(self._project_id, self._experiment_id, metrics=plan.metrics)
        except requests.RequestException as e:
            # 内容已上传但指标未确认：全部保留，重试时整流程重跑（重新签名 + 覆盖写）
            return UploadResult(
                remaining_records=list(records),
                failures=plan.blocked_failures
                + [transport_failure(plan.metric_records, e, source="swanlab_api", stage="submit_metrics")],
            )
        # 内容与指标均确认成功：完成文件行，确认产出过指标的 Record
        if self._tracker is not None:
            for payload in plan.payloads:
                self._tracker.finish_file(payload.tracker_key)
        return UploadResult(
            confirmed_records=plan.metric_records,
            remaining_records=[failure.records[0] for failure in plan.blocked_failures],
            failures=plan.blocked_failures,
        )

    # ── 文件保存上传 ──

    def upload_save(self, records: Sequence[Record]) -> UploadResult:
        """内部保存仅上传内容，用户保存（CUSTOM）走文件上传（小文件 presigned URL，大文件分片）。

        skip_store 下内部保存在 payload 中，CUSTOM 只读 source_path。
        """
        confirmed: List[Record] = []
        remaining: List[Record] = []
        failures: List[UploadFailure] = []
        custom_records: List[Record] = []
        for record in records:
            if not record.HasField("save"):
                failures.append(self._blocked_failure([record], "record without save payload"))
                remaining.append(record)
                continue
            save = record.save
            if save.type == SaveType.SAVE_TYPE_CUSTOM:
                # CUSTOM 的 payload 必须缺席，present 视为协议违约保留（防大文件误入内存通道）
                if save.HasField("payload"):
                    console.warning(f"CUSTOM save must not carry payload: {save.name or save.source_path}")
                    failures.append(self._blocked_failure([record], "CUSTOM save must not carry payload"))
                    remaining.append(record)
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
        """上传单条内部 save：读数据 → 解码 → 上传（payload 与磁盘双通道）。

        读取/解码失败与上传错误均转为结果内的失败，记录保留；不在此吞错。
        """
        save = record.save
        data = self.read_internal_save(save)
        if data is None:
            return UploadResult(
                remaining_records=[record],
                failures=[self._blocked_failure([record], f"internal save content unavailable: {save.source_path}")],
            )
        content = self.decode_internal_save(save.type, data)
        if content is None:
            return UploadResult(
                remaining_records=[record],
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
            return UploadResult(
                remaining_records=[record],
                failures=[transport_failure([record], e, source="swanlab_api", stage="submit_metrics")],
            )
        return UploadResult(confirmed_records=[record])

    def read_internal_save(self, save: SaveRecord) -> Optional[bytes]:
        """内部 save 数据源：payload 优先，默认模式回退磁盘；无法访问时告警并返回 None。"""
        if save.HasField("payload"):
            return save.payload
        if self._ctx.config.skip_store:
            console.warning(f"Internal save payload missing with skip_store enabled, skipping: {save.type}")
            return None
        source = resolve_save_source(self._ctx, save)
        if source is None:
            console.warning(f"Save file not found: {save.source_path}")
            return None
        try:
            return source.read_bytes()
        except OSError as e:
            console.warning(f"Failed to read save file, file kept at {source}: {e}")
            return None

    def decode_internal_save(self, save_type: SaveType, data: bytes) -> Optional[Any]:
        """按 SaveType 解码内部 save 字节：metadata/config 为 dict，requirements/conda 为非空 str。

        脏数据与空内容告警后返回 None，由调用方保留记录并报告失败。
        """
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
            console.warning(f"Failed to parse internal save payload: type={save_type}, error={e}")
            return None
        console.warning(f"Unknown internal save type {save_type}")
        return None


__all__ = [
    "HttpRecordSender",
    "encode_column",
]
