"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传抽象（HTTP 传输层）
"""

from __future__ import annotations

import json
import math
import mimetypes
from collections.abc import Callable, Sequence
from io import BytesIO
from pathlib import Path, PurePosixPath
from typing import Literal, Optional, Union, cast

import yaml

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.api.save import (
    complete_multipart_save,
    complete_save_files,
    prepare_multipart_save,
    prepare_save_files,
)
from swanlab.sdk.internal.core_python.api.upload import (
    upload_columns,
    upload_conda,
    upload_config,
    upload_log,
    upload_media,
    upload_metadata,
    upload_requirements,
    upload_resource,
    upload_saves,
    upload_scalar,
)
from swanlab.sdk.internal.core_python.context import CoreContext
from swanlab.sdk.internal.core_python.pkg.executor import SafeThreadPoolExecutor
from swanlab.sdk.internal.pkg import adapter, client, console, safe
from swanlab.sdk.internal.pkg.client.session import SessionWithRetry
from swanlab.sdk.typings.core_python.api.save import (
    CompletedMultipartSaveFile,
    CompletedPart,
    CompleteSaveFile,
    MultipartPart,
    PrepareMultipartSaveFile,
    PrepareSaveFile,
    SaveFileEntry,
)
from swanlab.sdk.typings.core_python.api.upload import (
    UploadColumn,
    UploadLog,
    UploadLogBatch,
    UploadMedia,
    UploadMediaBatch,
    UploadResource,
    UploadScalar,
    UploadScalarBatch,
)

from .helper import compute_md5


class HttpRecordSender:
    """HTTP 传输层上传 sender

    内部捕获 config、metadata、requirements、conda 等特殊文件的上传逻辑，
    并将其他类型的上传交给上层处理——这部分上传将有重试机制
    """

    def __init__(self, ctx: CoreContext, username: str, project: str, project_id: str, experiment_id: str) -> None:
        self._ctx = ctx
        self._username = username
        self._project = project
        self._project_id = project_id
        self._experiment_id = experiment_id
        # 资源上传 session，作用在Sender对象中复用TCP链接
        self._buffer_session: Optional[SessionWithRetry] = None
        self._upload_handlers: dict[str, Callable[[Sequence[Record]], None]] = {
            "column": self.upload_column,
            "scalar": self.upload_scalar,
            "media": self.upload_media,
            "config": self.upload_config,
            "log": self.upload_log,
            "metadata": self.upload_metadata,
            "requirements": self.upload_requirements,
            "conda": self.upload_conda,
            "save": self.upload_save,
        }

    def _ensure_session(self) -> SessionWithRetry:
        """懒初始化资源上传 session（复用 TCP 连接）。"""
        if self._buffer_session is None:
            self._buffer_session = client.session.create()
        return self._buffer_session

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
        paths: list[str] = []
        buffers: list[BytesIO] = []
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
                buffer_chunk: list[BytesIO] = []
                for media in media_record.value.items:
                    # 目前约定的本地文件路径格式为：media/<type>/filename，不区分key
                    # 约定保存到对象存储中的文件路径类似于本地文件路径
                    medium = adapter.medium[media_record.type]
                    remote_path = PurePosixPath("media", medium, media.filename)
                    local_path = self._ctx.media_dir / medium / media.filename
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
        upload_resource(self._ensure_session(), self._experiment_id, paths=paths, buffers=buffers)
        upload_media(self._project_id, self._experiment_id, metrics=metrics)

    def upload_log(self, records: Sequence[Record]) -> None:
        with safe.block(message="Failed to upload terminal logs, skipping"):
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
            upload_log(self._project_id, self._experiment_id, metrics=metrics)

    def upload_config(self, _: Sequence[Record]) -> None:
        with safe.block(message=f"Failed to upload config, skipping; file kept at {self._ctx.config_file}"):
            with open(self._ctx.config_file, "r", encoding="utf-8") as f:
                content = yaml.safe_load(f)
            if isinstance(content, dict):
                upload_config(self._username, self._project, self._experiment_id, content=content)

    def upload_metadata(self, _: Sequence[Record]) -> None:
        with safe.block(message=f"Failed to upload metadata, skipping; file kept at {self._ctx.metadata_file}"):
            with open(self._ctx.metadata_file, "r", encoding="utf-8") as f:
                content = json.load(f)
            if isinstance(content, dict):
                upload_metadata(self._username, self._project, self._experiment_id, content=content)

    def upload_requirements(self, _: Sequence[Record]) -> None:
        with safe.block(message=f"Failed to upload requirements, skipping; file kept at {self._ctx.requirements_file}"):
            with open(self._ctx.requirements_file, "r", encoding="utf-8") as f:
                content = f.read()
            if len(content) > 0:
                upload_requirements(self._username, self._project, self._experiment_id, content=content)

    def upload_conda(self, _: Sequence[Record]) -> None:
        with safe.block(message=f"Failed to upload conda, skipping; file kept at {self._ctx.conda_file}"):
            with open(self._ctx.conda_file, "r", encoding="utf-8") as f:
                content = f.read()
            if len(content) > 0:
                upload_conda(self._username, self._project, self._experiment_id, content=content)

    # ── 文件保存上传 ──

    def upload_save(self, records: Sequence[Record]) -> None:
        """处理 SaveRecord 上传：小文件走 presigned URL，大文件走分片上传。"""
        config = self._ctx.config
        # 1. 收集合法文件（不读内容，避免大批次内存爆炸）
        pending: list[tuple[str, int, str]] = []  # (source_path, size, name)
        for record in records:
            if not record.HasField("save"):
                continue
            save = record.save
            source = Path(save.source_path)
            if not source.is_file():
                console.warning(f"Save file not found, skipping: {save.source_path}")
                continue
            size: Optional[int] = None
            with safe.block(message=f"Failed to stat save file, skipping: {save.source_path}"):
                size = source.stat().st_size
            if size is None:
                continue
            if size > config.save_size:
                console.warning(
                    f"Save file exceeds size limit ({size} > {config.save_size}), skipping: {save.source_path}"
                )
                continue
            pending.append((save.source_path, size, save.name))

        if not pending:
            console.warning("No valid files to save.")
            return

        # 2. 校验批次限制：数量超限时直接拒绝整个批次
        # FIXME: pending过长时根据save_batch分片，而非直接拒绝
        if len(pending) > config.save_batch:
            console.warning(
                f"Save batch size ({len(pending)}) exceeds limit ({config.save_batch}), skipping entire batch"
            )
            return

        # 3. 异步计算 MD5 + 按大小分流上传
        with safe.block(message="Failed to upload save files, skipping"):
            file_entries: list[SaveFileEntry] = []
            with SafeThreadPoolExecutor(max_workers=4) as executor:
                md5_futures = {source_path: executor.submit(compute_md5, source_path) for source_path, _, _ in pending}
                for source_path, size, name in pending:
                    md5 = md5_futures[source_path].result()
                    if not isinstance(md5, str):
                        console.warning(f"Save file MD5 computation failed, skipping: {source_path}")
                        continue
                    mime_type = mimetypes.guess_type(source_path)[0] or "application/octet-stream"
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
                self._upload_small_saves(small_files)
            if large_files:
                self._upload_large_saves(large_files)

    def _upload_small_saves(self, files: list[SaveFileEntry]) -> None:
        """小文件：prepare → presigned PUT（并发）→ complete（逐文件状态）。"""
        prepare_items: list[PrepareSaveFile] = [
            {"path": f["path"], "size": f["size"], "md5": f["md5"], "mimeType": f["mime_type"]} for f in files
        ]
        resp = prepare_save_files(self._experiment_id, files=prepare_items)
        urls = resp["urls"]
        if len(urls) != len(files):
            console.warning(f"Save prepare returned {len(urls)} URLs for {len(files)} files, skipping")
            return
        session = self._ensure_session()
        resources: list[UploadResource] = [
            {"url": url, "source_path": f["source_path"], "content_type": f["mime_type"]} for f, url in zip(files, urls)
        ]
        source_to_path = {f["source_path"]: f["path"] for f in files}
        statuses: dict[str, str] = {f["path"]: "FAILED" for f in files}
        with safe.block(message="Failed to upload save files, skipping"):
            uploaded_paths = upload_saves(session, resources=resources)
            for source_path in uploaded_paths:
                statuses[source_to_path[source_path]] = "UPLOADED"

        complete_items: list[CompleteSaveFile] = [
            {"path": path, "state": cast(Literal["UPLOADED", "FAILED"], state)} for path, state in statuses.items()
        ]
        with safe.block(message="Failed to complete save upload"):
            complete_save_files(self._experiment_id, files=complete_items)
        synced = sum(1 for s in statuses.values() if s == "UPLOADED")
        if synced > 0:
            console.info(f"Synced {synced} save file(s)")

    def _upload_large_saves(self, files: list[SaveFileEntry]) -> None:
        """大文件：分片上传。"""
        part_size = self._ctx.config.save_part
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
        resp = prepare_multipart_save(self._experiment_id, files=prepare_items)
        multipart_files = resp["files"]
        if len(multipart_files) != len(files):
            console.warning(
                f"Multipart prepare returned {len(multipart_files)} entries for {len(files)} files, skipping"
            )
            return
        completed: list[CompletedMultipartSaveFile] = []
        failed: list[CompleteSaveFile] = []
        for file_info, multipart_info in zip(files, multipart_files):
            upload_id = multipart_info["uploadId"]
            parts = multipart_info["parts"]
            self._ensure_session()
            completed_parts = self._upload_multipart_parts(file_info, parts)
            if completed_parts is None:
                console.warning(f"Multipart upload failed for {file_info['path']}, skipping")
                failed.append({"path": file_info["path"], "state": "FAILED"})
                continue
            completed.append(
                {
                    "path": file_info["path"],
                    "parts": completed_parts,
                    "uploadId": upload_id,
                }
            )
        if completed:
            with safe.block(message="Failed to complete multipart save upload"):
                complete_multipart_save(self._experiment_id, files=completed)
            console.info(f"Synced {len(completed)} save file(s)")
        if failed:
            with safe.block(message="Failed to report multipart save failures"):
                complete_save_files(self._experiment_id, files=failed)

    def _upload_multipart_parts(
        self, file_info: SaveFileEntry, parts: Sequence[MultipartPart]
    ) -> Optional[list[CompletedPart]]:
        """上传单个大文件的所有分片，返回 CompletedPart 列表。"""
        assert self._buffer_session is not None
        completed_parts: list[CompletedPart] = []
        source_path = file_info["source_path"]
        part_size = self._ctx.config.save_part
        session = self._buffer_session

        def upload_part(part_info: MultipartPart) -> CompletedPart:
            part_number = part_info["partNumber"]
            url = part_info["url"]
            offset = (part_number - 1) * part_size
            with open(source_path, "rb") as f:
                f.seek(offset)
                data = f.read(part_size)
            resp = session.put(url, data=data, headers={"Content-Type": "application/octet-stream"})
            resp.raise_for_status()
            etag = resp.headers.get("ETag", "").strip('"')
            return {"partNumber": part_number, "etag": etag}

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
