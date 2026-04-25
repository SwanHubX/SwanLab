"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传抽象（HTTP 传输层）
"""

from __future__ import annotations

import json
from collections.abc import Callable, Sequence
from pathlib import Path

import yaml

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.api.upload import (
    upload_conda,
    upload_config,
    upload_console,
    upload_metadata,
    upload_requirements,
)
from swanlab.sdk.internal.pkg import adapter, console, safe
from swanlab.sdk.typings.core_python.api.upload import UploadLog, UploadLogMetrics


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
        self._console_epoch = 1

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
        handler(records)

    def upload_column(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_column (request mapping pending).")

    def upload_scalar(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_scalar (request mapping pending).")

    def upload_media(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_media (request mapping pending).")

    def upload_console(self, records: Sequence[Record]) -> None:
        console_records = [r.console for r in records if r.HasField("console")]
        if not console_records:
            return
        with safe.block(message="Failed to upload console logs, skipping"):
            metrics: UploadLogMetrics = []
            for console_record in console_records:
                if console_record.HasField("timestamp"):
                    create_time = console_record.timestamp.ToJsonString()
                    metric: UploadLog = {
                        "level": adapter.level[console_record.stream],
                        "epoch": self._console_epoch,
                        "message": console_record.line,
                        "create_time": create_time,
                    }
                    metrics.append(metric)
                self._console_epoch += 1
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


__all__ = [
    "HttpRecordSender",
]
