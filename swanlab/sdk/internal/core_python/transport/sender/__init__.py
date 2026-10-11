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
from typing import TYPE_CHECKING, Any, Literal, Optional, Union, cast

import yaml

from swanlab.exceptions import ApiError
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

from .files import put_media_resources, upload_save_files
from .utils import collect_media_payloads, encode_column, resolve_save_source

if TYPE_CHECKING:
    from swanlab.sdk.internal.core_python.transport.tracker import UploadTracker
    from swanlab.sdk.internal.pkg.client.session import SessionWithRetry


class HttpRecordSender:
    """Record 与后端对接的真实上传接口。

    全部 upload_* 方法显式定义于本类，每个方法是一个完整短流程：
    编码/收集（utils.py）→ 对象存储字节传输（files.py）→ 业务指标提交（api 层）。
    本类只持有三个实例状态：ctx、懒初始化的上传 session、tracker 及其去重集合；
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
        self._upload_handlers: dict[str, Callable[[Sequence[Record]], None]] = {
            "column": self.upload_column,
            "scalar": self.upload_scalar,
            "media": self.upload_media,
            "log": self.upload_log,
            "save": self.upload_save,
        }

    def set_tracker(self, tracker: UploadTracker) -> None:
        self._tracker = tracker

    def _ensure_session(self) -> SessionWithRetry:
        """懒初始化资源上传 session（复用 TCP 连接）。"""
        if self._buffer_session is None:
            self._buffer_session = client.session.create()
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

    # ── 通用上传入口 ──

    def upload(self, record_type: str, records: Sequence[Record]) -> None:
        """通用上传入口。"""
        if not records:
            return
        handler = self._upload_handlers.get(record_type)
        if handler is None:
            console.warning(
                f"No upload handler registered for record type {record_type!r}; skipping {len(records)} record(s)."
            )
            return
        try:
            handler(records)
        except ApiError as e:
            if e.response.status_code >= 500:
                # 基础设施错误（断网/5xx），向上抛出交由 transport 重试；raise 路径不会执行 advance
                raise
            # 4xx：后端业务拒绝，重试也不会成功；不计入 uploaded，进度条原地不动以如实反映失败
            console.warning(f"Failed to upload {record_type} records, skipping; error: {e}")
            return
        # 仅当本批真正上传成功后才递进 uploaded（网络失败会 raise，走不到这里）
        if self._tracker is not None:
            self._tracker.advance_records(len(records))

    # ── 普通记录 ──

    def upload_column(self, records: Sequence[Record], batch_size: int = 3000) -> None:
        columns = [encode_column(record.column) for record in records]
        if not columns:
            return
        for i in range(0, len(columns), batch_size):
            batch = columns[i : i + batch_size]
            upload_columns(self._username, self._project, columns={"series": batch})

    def upload_scalar(self, records: Sequence[Record]) -> None:
        metrics: UploadScalarBatch = []
        for record in records:
            if not record.HasField("scalar"):
                continue
            scalar_record = record.scalar
            if scalar_record.HasField("timestamp"):
                create_time = scalar_record.timestamp.ToJsonString()
                data = scalar_record.value.number if math.isfinite(scalar_record.value.number) else "nan"
                metric: UploadScalar = {
                    "key": scalar_record.key,
                    "index": scalar_record.step,
                    "data": cast(Union[float, Literal["nan"]], data),
                    "create_time": create_time,
                }
                metrics.append(metric)
        console.debug(f"HTTP upload: upload_column with metrics count: {len(metrics)}")
        upload_scalar(self._project_id, self._experiment_id, metrics=metrics)

    def upload_log(self, records: Sequence[Record]) -> None:
        console.debug(f"Preparing log HTTP request: records={len(records)}", write_to_tty=False)
        metrics: UploadLogBatch = []
        for record in records:
            if not record.HasField("log"):
                continue
            log_record = record.log
            if log_record.HasField("timestamp"):
                create_time = log_record.timestamp.ToJsonString()
                metric: UploadLog = {
                    "level": adapter.level[log_record.level],
                    "epoch": log_record.epoch,
                    "message": log_record.line,
                    "create_time": create_time,
                }
                metrics.append(metric)
        console.debug(f"Sending log HTTP request: metrics={len(metrics)}", write_to_tty=False)
        upload_log(self._project_id, self._experiment_id, metrics=metrics)
        console.debug(f"Log HTTP request completed: metrics={len(metrics)}", write_to_tty=False)

    # ── 媒体 ──

    def upload_media(self, records: Sequence[Record]) -> None:
        """媒体：收集载荷 → 对象存储 PUT → 提交媒体指标，两阶段均成功才算完成。"""
        plan = collect_media_payloads(self._ctx, records)
        if not plan.payloads:
            return
        put_media_resources(self._ensure_session(), self._experiment_id, plan.payloads, tracker=self._tracker)
        upload_media(self._project_id, self._experiment_id, metrics=plan.metrics)

    # ── 文件保存上传 ──

    def upload_save(self, records: Sequence[Record]) -> None:
        """内部保存仅上传内容，用户保存（CUSTOM）走文件上传（小文件 presigned URL，大文件分片）。

        skip_store 下内部保存在 payload 中，CUSTOM 只读 source_path。
        """
        # 1. 分流：内部保存直接上传，CUSTOM 收集后走文件上传
        save_records: list[Record] = []
        for record in records:
            if not record.HasField("save"):
                continue
            save = record.save
            if save.type == SaveType.SAVE_TYPE_CUSTOM:
                # CUSTOM 的 payload 必须缺席，present 视为协议违约丢弃（防大文件误入内存通道）
                if save.HasField("payload"):
                    console.warning(f"CUSTOM save must not carry payload, skipping: {save.name or save.source_path}")
                    continue
                save_records.append(record)
            else:
                self.upload_internal_save(record)
        if not save_records:
            console.debug("No user save records to upload after filtering; all saves are internal files.")
            return
        # 2. CUSTOM 走对象存储文件上传流程
        upload_save_files(
            self._ctx,
            self._ensure_session(),
            save_records,
            tracker=self._tracker,
            on_part_progress=self._mark_upload_progress,
        )

    def upload_internal_save(self, record: Record) -> None:
        """上传单条内部 save：读数据 → 解码 → 上传（payload 与磁盘双通道）。

        前两步返回 None 即已告警跳过；上传异常不捕获，5xx 交 Transport 重试、4xx 由外层跳过。
        """
        save = record.save
        data = self.read_internal_save(save)
        if data is None:
            return
        content = self.decode_internal_save(save.type, data)
        if content is None:
            return
        if save.type == SaveType.SAVE_TYPE_METADATA:
            upload_metadata(self._username, self._project, self._experiment_id, content=content)
        elif save.type == SaveType.SAVE_TYPE_REQUIREMENTS:
            upload_requirements(self._username, self._project, self._experiment_id, content=content)
        elif save.type == SaveType.SAVE_TYPE_CONDA:
            upload_conda(self._username, self._project, self._experiment_id, content=content)
        elif save.type == SaveType.SAVE_TYPE_CONFIG:
            upload_config(self._username, self._project, self._experiment_id, content=content)

    def read_internal_save(self, save: SaveRecord) -> Optional[bytes]:
        """内部 save 数据源：payload 优先，默认模式回退磁盘；无法访问时告警并返回 None。"""
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
        """按 SaveType 解码内部 save 字节：metadata/config 为 dict，requirements/conda 为非空 str。

        脏数据告警后返回 None，空内容静默返回 None，均不进入 Transport 无上限重试。
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
            console.warning(f"Failed to parse internal save payload, skipping: type={save_type}, error={e}")
            return None
        console.warning(f"Unknown internal save type {save_type}, skipping")
        return None


__all__ = [
    "HttpRecordSender",
    "encode_column",
]
