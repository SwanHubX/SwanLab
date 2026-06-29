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
import threading
from collections.abc import Callable, Sequence
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import IO, TYPE_CHECKING, Literal, Optional, Union, cast

import yaml

from swanlab.exceptions import ApiError
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnClass, ColumnRecord, ColumnType
from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
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
from swanlab.sdk.internal.core_python.utils import ProgressFileWrapper, get_buffer_size
from swanlab.sdk.internal.pkg import adapter, client, console, safe
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

from .helper import compute_md5, save_tracker_key

if TYPE_CHECKING:
    from .tracker import UploadTracker


class HttpRecordSender:
    """HTTP 传输层上传 sender

    内部捕获 config、metadata、requirements、conda 等特殊文件的上传逻辑，
    并将其他类型的上传交给上层处理——这部分上传将有重试机制
    """

    def __init__(self, ctx: CoreContext) -> None:
        self._ctx = ctx
        self._username = ctx.username
        self._project = ctx.project
        self._project_id = ctx.project_id
        self._experiment_id = ctx.experiment_id
        self._tracker: Optional[UploadTracker] = None
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

    def _track_file(self, key: str, path: str, size: int) -> None:
        """在 tracker 中注册文件显示行，重复 key 由 UploadTracker 自动去重。"""
        if self._tracker is None or size <= 0:
            return
        self._tracker.add_file(key=key, path=path, total=size)

    def _mark_upload_progress(self, file_key: str, upload_key: str, current_bytes: int, *, final: bool = False) -> None:
        """推送字节进度到 tracker，跳过已完成的上传事件 key。

        final=True 时将 upload_key 标记为已完成，后续对该 key 的调用（包括中间进度）将被跳过。
        """
        with self._completed_upload_keys_lock:
            if upload_key in self._completed_upload_keys:
                return
            if final:
                self._completed_upload_keys.add(upload_key)
        if self._tracker is not None:
            self._tracker.update_file_progress(file_key, upload_key, current_bytes)

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
                # 基础设施错误（断网/5xx），向上抛出交由 transport 重试；raise 路径不会执行 advance
                raise
            # 4xx：后端业务拒绝，重试也不会成功；不计入 uploaded，进度条原地不动以如实反映失败
            console.warning(f"Failed to upload {record_type} records, skipping; error: {e}")
            return
        # 仅当本批真正上传成功后才递进 uploaded（网络失败会 raise，走不到这里）
        if self._tracker is not None:
            self._tracker.advance_records(len(records))

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
        buffers: list[Union[IO[bytes], str, Path]] = []
        content_types: list[str] = []
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
                record_paths = []
                for media in media_record.value.items:
                    # 目前约定的本地文件路径格式为：media/<type>/filename，不区分key
                    # 约定保存到对象存储中的文件路径类似于本地文件路径
                    medium = adapter.medium[media_record.type]
                    remote_path = PurePosixPath("media", medium, media.filename)
                    local_path = self._ctx.media_dir / medium / media.filename

                    if not local_path.is_file():
                        continue

                    with safe.block(message="Failed to process media file, skipping"):
                        size = get_buffer_size(local_path)
                        remote_path_str = remote_path.as_posix()
                        tracker_key = f"{remote_path_str}:{size}"
                        self._track_file(tracker_key, local_path.as_posix(), size)
                        # 先计算 mime_type, 确保后续 record_paths/paths/buffers/content_types 原子化追加,
                        # 避免 safe.block 内中途异常导致列表长度不一致、MIME 类型错位
                        mime_type = mimetypes.guess_type(str(local_path))[0] or "application/octet-stream"
                        record_paths.append(remote_path_str)
                        paths.append(remote_path_str)
                        buffers.append(local_path)
                        content_types.append(mime_type)
                        if media.caption:
                            metric_chunk["more"].append({"caption": media.caption})
                if len(record_paths) > 0:
                    metric_chunk["data"] = record_paths
                    metrics.append(metric_chunk)
        if not (len(paths) == len(buffers) == len(content_types)):
            console.warning(
                f"Failed to upload media, quantity mismatch, skipping; "
                f"got {len(paths)} paths, {len(buffers)} buffers, {len(content_types)} content_types"
            )
            return
        # 先上传资源再上传媒体
        upload_resource(
            self._ensure_session(),
            self._experiment_id,
            paths=paths,
            buffers=buffers,
            content_types=content_types,
            tracker=self._tracker,
        )

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

    # ── 文件保存上传 ──

    def _resolve_save_source(self, save: SaveRecord) -> Optional[Path]:
        """解析 save 记录对应的可读本地文件路径。

        ``source_path`` 为训练机绝对路径；online/local/offline 本机运行时可读。
        sync 在另一台机器、以不同挂载根读取同一 run 目录时该绝对路径不可读，
        此时回退到当前 run 目录的 ``files`` 子目录按文件名重新定位：

        - 用户保存（CUSTOM）：按记录相对名 ``save.name`` 定位（镜像位置，与 ``create_save_links`` 约定一致）；
        - 内部保存（metadata/requirements/conda/config）：由 probe 直接写入 files 目录（真实文件），
          按 ``source_path`` 的 basename 定位。注意 config 的 ``save.name`` 为 ``"config"`` 而非 ``config.yaml``，
          故此处一律用 basename 而非 name。

        :param save: 文件保存记录
        :return: 可读的本地文件路径；若原始路径与回退路径均不可读则返回 None
        """
        primary = Path(save.source_path)
        if primary.is_file():
            return primary
        # source_path/name 可能由异构系统写入（例如训练在 Windows，sync 在 POSIX），
        # 其分隔符为反斜杠时，POSIX 的 Path 无法正确切分。用 PureWindowsPath 仅做分隔符解析，
        # 再交给本地 Path 复原，确保 basename 取正确文件名、子目录层级被正确还原。
        if save.type == SaveType.SAVE_TYPE_CUSTOM:
            fallback = self._ctx.files_dir / Path(*PureWindowsPath(save.name).parts)
        else:
            fallback = self._ctx.files_dir / PureWindowsPath(save.source_path).name
        if fallback.is_file():
            console.debug(f"Save source path not readable, recovered from run directory: {fallback}")
            return fallback
        return None

    def upload_save(self, records: Sequence[Record]) -> None:
        """
        内部保存（如config、metadata等）和用户保存（文件保存）共用 SaveRecord 结构
        目前它们在产品设计上暂未统一，换句话说上传config、metadata文件的时候并不会保存对应文件
        因此需要区分两者，内部保存仅上传，用户保存走save逻辑：小文件走 presigned URL，大文件走分片上传。
        """
        # 1. 区分内部保存和用户保存，过滤出需要走文件上传逻辑的记录，并处理内部保存的上传
        save_records = []
        for record in records:
            if not record.HasField("save"):
                continue
            save = record.save
            # 根据约定的 type 字段区分内部保存和用户保存，内部保存直接上传内容，用户保存走后续文件上传逻辑
            if save.type == SaveType.SAVE_TYPE_CUSTOM:
                save_records.append(record)
                continue
            # 内部保存（metadata/requirements/conda/config）：source_path 为训练机绝对路径，
            # 跨挂载根 sync 时不可读，回退到 run 目录 files 子目录按 basename 重新定位
            source_ref_path = self._resolve_save_source(save)
            if source_ref_path is None:
                console.warning(f"Save file not found, skipping: {save.source_path}")
                continue
            if save.type == SaveType.SAVE_TYPE_METADATA:
                with safe.block(message=f"Failed to upload metadata, skipping; file kept at {source_ref_path}"):
                    with open(source_ref_path, "r", encoding="utf-8") as f:
                        content = json.load(f)
                    if isinstance(content, dict):
                        upload_metadata(self._username, self._project, self._experiment_id, content=content)
            elif save.type == SaveType.SAVE_TYPE_REQUIREMENTS:
                with safe.block(message=f"Failed to upload requirements, skipping; file kept at {source_ref_path}"):
                    with open(source_ref_path, "r", encoding="utf-8") as f:
                        content = f.read()
                    if len(content) > 0:
                        upload_requirements(self._username, self._project, self._experiment_id, content=content)
            elif save.type == SaveType.SAVE_TYPE_CONDA:
                with safe.block(message=f"Failed to upload conda, skipping; file kept at {source_ref_path}"):
                    with open(source_ref_path, "r", encoding="utf-8") as f:
                        content = f.read()
                    if len(content) > 0:
                        upload_conda(self._username, self._project, self._experiment_id, content=content)
            elif save.type == SaveType.SAVE_TYPE_CONFIG:
                with safe.block(message=f"Failed to upload config, skipping; file kept at {source_ref_path}"):
                    with open(source_ref_path, "r", encoding="utf-8") as f:
                        content = yaml.safe_load(f)
                    if isinstance(content, dict):
                        upload_config(self._username, self._project, self._experiment_id, content=content)
            else:
                console.warning(f"Unknown save type {save.type} for record num {record.num}, skipping")
        if not save_records:
            console.debug("No user save records to upload after filtering; all saves are internal files.")
            return

        # 2. save 逻辑
        config = self._ctx.config
        # 2.1 收集合法文件（不读内容，避免大批次内存爆炸）
        pending: list[tuple[str, int, str]] = []  # (source_path, size, name)
        for record in save_records:
            if not record.HasField("save"):
                continue
            save = record.save
            source_ref_path = self._resolve_save_source(save)
            if source_ref_path is None:
                console.warning(f"Save file not found, skipping: {save.source_path}")
                continue
            size: Optional[int] = None
            with safe.block(message=f"Failed to stat save file, skipping: {source_ref_path}"):
                size = source_ref_path.stat().st_size
            if size is None:
                continue
            if size > config.save_size:
                console.warning(
                    f"Save file exceeds size limit ({size} > {config.save_size}), skipping: {source_ref_path}"
                )
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

        # 2.2 按 save_batch 拆分文件列表，避免单次 prepare/complete 超过后端限制
        for index in range(0, len(pending), config.save_batch):
            self._upload_save_batch(pending[index : index + config.save_batch])

    def _upload_save_batch(self, pending: Sequence[tuple[str, int, str]]) -> None:
        """上传一个 save_batch 内的文件列表。"""
        config = self._ctx.config
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
        resources: list[UploadResource] = []
        for f, url in zip(files, urls):
            tracker_key = save_tracker_key(f)
            self._track_file(tracker_key, f["source_path"], f["size"])
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
                tracker=self._tracker,
            )
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
            file_key = save_tracker_key(file_info)
            self._track_file(file_key, file_info["source_path"], file_info["size"])
            self._ensure_session()
            completed_parts = self._upload_multipart_parts(file_info, parts)
            if completed_parts is None:
                console.warning(f"Multipart upload failed for {file_info['path']}, skipping")
                failed.append({"path": file_info["path"], "state": "FAILED"})
                continue
            if self._tracker is not None:
                self._tracker.finish_file(file_key)
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
        self,
        file_info: SaveFileEntry,
        parts: Sequence[MultipartPart],
    ) -> Optional[list[CompletedPart]]:
        """上传单个大文件的所有分片，返回 CompletedPart 列表。"""
        assert self._buffer_session is not None
        completed_parts: list[CompletedPart] = []
        source_path = file_info["source_path"]
        part_size = self._ctx.config.save_part
        session = self._buffer_session
        file_key = save_tracker_key(file_info)

        def upload_part(part_info: MultipartPart) -> CompletedPart:
            part_number = part_info["partNumber"]
            url = part_info["url"]
            offset = (part_number - 1) * part_size
            part_key = f"{file_key}:part:{part_number}"
            part_length = min(part_size, max(file_info["size"] - offset, 0))

            def on_read(current_bytes):
                self._mark_upload_progress(file_key, part_key, current_bytes)

            try:
                with open(source_path, "rb") as f:
                    wrapped_data = ProgressFileWrapper(f, on_read, part_length, offset=offset, size=part_length)
                    resp = session.put(url, data=wrapped_data, headers={"Content-Type": "application/octet-stream"})
                resp.raise_for_status()
                self._mark_upload_progress(file_key, part_key, part_length, final=True)
                etag = resp.headers.get("ETag", "").strip('"')
                return {"partNumber": part_number, "etag": etag}
            except Exception:
                self._mark_upload_progress(file_key, part_key, 0)
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
            y_range = record.y_range
            min_value = y_range.min if y_range.HasField("min") else None
            max_value = y_range.max if y_range.HasField("max") else None
            column["yRange"] = (min_value, max_value)
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
