"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传抽象（HTTP 传输层）
"""

from __future__ import annotations

import json
import math
from _io import BytesIO
from collections.abc import Callable, Sequence
from pathlib import Path, PurePosixPath
from typing import List, Literal, Optional, Union, cast

import yaml
from requests.sessions import Session

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.api.upload import (
    upload_columns,
    upload_conda,
    upload_config,
    upload_console,
    upload_media,
    upload_metadata,
    upload_requirements,
    upload_resource,
    upload_scalar,
)
from swanlab.sdk.internal.pkg import adapter, client, console, safe
from swanlab.sdk.typings.core_python.api.upload import (
    UploadColumn,
    UploadLog,
    UploadLogBatch,
    UploadMedia,
    UploadMediaBatch,
    UploadScalar,
    UploadScalarBatch,
)


class HttpRecordSender:
    """HTTP 传输层上传 sender

    内部捕获 config、metadata、requirements、conda 等特殊文件的上传逻辑，
    并将其他类型的上传交给上层处理——这部分上传将有重试机制
    """

    def __init__(self, run_dir: Path, username: str, project: str, project_id: str, experiment_id: str) -> None:
        self._run_dir = run_dir
        self._username = username
        self._project = project
        self._project_id = project_id
        self._experiment_id = experiment_id
        # 资源上传 session，作用在Sender对象中复用TCP链接
        self._buffer_session: Optional[Session] = None
        self._upload_handlers: dict[str, Callable[[Sequence[Record]], None]] = {
            "column": self.upload_column,
            "scalar": self.upload_scalar,
            "media": self.upload_media,
            "config": self.upload_config,
            "console": self.upload_console,
            "metadata": self.upload_metadata,
            "requirements": self.upload_requirements,
            "conda": self.upload_conda,
        }

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
                # 基础设施错误，不按业务逻辑处理，此时直接抛出异常，交给上游重试
                raise
            console.warning(f"Failed to upload {record_type} records, skipping; error: {e}")

    def upload_column(self, records: Sequence[Record], batch_size: int = 3000) -> None:
        columns = []
        for record in records:
            r = encode_column(record.column)
            if r:
                columns.append(r)
        if not columns:
            return
        for i in range(0, len(columns), batch_size):
            upload_columns(self._experiment_id, columns=columns[i : i + batch_size])

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

    def upload_media(self, records: Sequence[Record]) -> None:
        metrics: UploadMediaBatch = []
        paths: List[str] = []
        buffers: List[BytesIO] = []
        for record in records:
            if not record.HasField("media"):
                continue
            media_record = record.media
            if media_record.HasField("timestamp"):
                create_time = media_record.timestamp.ToJsonString()
                metric_chunk: UploadMedia = {
                    "key": media_record.key,
                    "index": media_record.step,
                    "data": [],
                    "more": [],
                    "create_time": create_time,
                }
                buffer_chunk: List[BytesIO] = []
                for media in media_record.value.items:
                    # 目前约定的本地文件路径格式为：media/<type>/filename，不区分key
                    # 约定保存到对象存储中的文件路径类似
                    remote_path = PurePosixPath("media", adapter.medium[media_record.type], media.filename)
                    local_path = self._run_dir.joinpath(*remote_path.parts)
                    # 读取本地文件内容
                    buffer: Optional[BytesIO] = None
                    with safe.block(message=f"Failed to read local file: {local_path}, skipping"):
                        with open(local_path, "rb") as f:
                            buffer = BytesIO(f.read())
                    if buffer is None:
                        continue
                    buffer_chunk.append(buffer)
                    metric_chunk["data"].append(remote_path.as_posix())
                    if media.caption:
                        metric_chunk["more"].append({"caption": media.caption})
                if len(metric_chunk["data"]) > 0:
                    metrics.append(metric_chunk)
                    paths.extend(metric_chunk["data"])
                    buffers.extend(buffer_chunk)
        if len(paths) != len(buffers):
            console.warning(
                f"Failed to upload media, quantity mismatch, skipping; got {len(paths)} paths, {len(buffers)} buffers"
            )
            return
        # 先上传资源再上传媒体
        if not self._buffer_session:
            self._buffer_session = client.session.create(default_retry=2)
        upload_resource(self._buffer_session, self._experiment_id, paths=paths, buffers=buffers)
        upload_media(self._project_id, self._experiment_id, metrics=metrics)

    def upload_console(self, records: Sequence[Record]) -> None:
        with safe.block(message="Failed to upload console logs, skipping"):
            metrics: UploadLogBatch = []
            for record in records:
                if not record.HasField("console"):
                    continue
                console_record = record.console
                if console_record.HasField("timestamp"):
                    create_time = console_record.timestamp.ToJsonString()
                    metric: UploadLog = {
                        "level": adapter.level[console_record.stream],
                        "epoch": console_record.epoch,
                        "message": console_record.line,
                        "create_time": create_time,
                    }
                    metrics.append(metric)
            upload_console(self._project_id, self._experiment_id, metrics=metrics)

    def upload_config(self, _: Sequence[Record]) -> None:
        config_path = self._run_dir / "files" / "config.yaml"
        with safe.block(message=f"Failed to upload config, skipping; file kept at {config_path}"):
            with open(config_path, "r", encoding="utf-8") as f:
                content = yaml.safe_load(f)
            if isinstance(content, dict):
                upload_config(self._username, self._project, self._experiment_id, content=content)

    def upload_metadata(self, _: Sequence[Record]) -> None:
        metadata_path = self._run_dir / "files" / "swanlab-metadata.json"
        with safe.block(message=f"Failed to upload metadata, skipping; file kept at {metadata_path}"):
            with open(metadata_path, "r", encoding="utf-8") as f:
                content = json.load(f)
            if isinstance(content, dict):
                upload_metadata(self._username, self._project, self._experiment_id, content=content)

    def upload_requirements(self, _: Sequence[Record]) -> None:
        requirements_path = self._run_dir / "files" / "requirements.txt"
        with safe.block(message=f"Failed to upload requirements, skipping; file kept at {requirements_path}"):
            with open(requirements_path, "r", encoding="utf-8") as f:
                content = f.read()
            if len(content) > 0:
                upload_requirements(self._username, self._project, self._experiment_id, content=content)

    def upload_conda(self, _: Sequence[Record]) -> None:
        conda_path = self._run_dir / "files" / "conda.yaml"
        with safe.block(message=f"Failed to upload conda, skipping; file kept at {conda_path}"):
            with open(conda_path, "r", encoding="utf-8") as f:
                content = f.read()
            if len(content) > 0:
                upload_conda(self._username, self._project, self._experiment_id, content=content)


def encode_column(record: ColumnRecord) -> Optional[UploadColumn]:
    """
    将列记录编码为后端所需的格式（DTO）
    """
    column: UploadColumn = {"key": record.column_key, "type": adapter.column[record.column_type]}
    # class: 是否为系统列
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["class"] = cast(Literal["SYSTEM"], "SYSTEM")
    # name: 列的展示名称
    if record.column_name:
        column["name"] = record.column_name
    # section name: section 名称
    if record.section_name:
        column["sectionName"] = record.section_name
    # section type: section 类型，目前仅处理系统列
    if record.column_class == ColumnClass.COLUMN_CLASS_SYSTEM:
        column["sectionType"] = cast(Literal["SYSTEM"], "SYSTEM")
    # yRange: 数值列的 y 轴范围
    if record.column_type == ColumnType.COLUMN_TYPE_SCALAR:
        if record.HasField("y_range"):
            column["yRange"] = (record.y_range.min, record.y_range.max)
    # chartName: 图表名称
    if record.chart_name:
        column["chartName"] = record.chart_name
    # chartIndex: 图表索引
    if record.chart_index:
        column["chartIndex"] = record.chart_index
    # metricName: 指标名称
    if record.metric_name:
        column["metricName"] = record.metric_name
    # metricColors: 指标颜色
    if record.HasField("metric_colors"):
        column["metricColors"] = (record.metric_colors.light, record.metric_colors.dark)
    return column


__all__ = [
    "HttpRecordSender",
]
