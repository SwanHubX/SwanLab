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
from swanlab.sdk.internal.core_python.api.environment import (
    upload_conda,
    upload_config,
    upload_metadata,
    upload_requirements,
)
from swanlab.sdk.internal.pkg import console, safe


class HttpRecordSender:
    """HTTP 传输层上传 sender

    内部捕获 config、metadata、requirements、conda 等特殊文件的上传逻辑，
    并将其他类型的上传交给上层处理——这部分上传将有重试机制
    """

    def __init__(self, run_dir: Path, username: str, project: str, cuid: str) -> None:
        self._run_dir = run_dir
        self._username = username
        self._project = project
        self._cuid = cuid
        self._upload_handlers: dict[str, Callable[[Sequence[Record]], None]] = {
            "column": self.upload_column,
            "data": self.upload_data,
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
        handler(records)

    def upload_column(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_column (request mapping pending).")

    def upload_data(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_data (request mapping pending).")

    def upload_console(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_console (request mapping pending).")

    def upload_config(self, _: Sequence[Record]) -> None:
        config_path = self._run_dir / "files" / "config.yaml"
        with safe.block(message=f"Failed to upload config, skipping; file kept at {config_path}"):
            with open(config_path, "r", encoding="utf-8") as f:
                content = yaml.safe_load(f)
            if isinstance(content, dict):
                upload_config(self._username, self._project, self._cuid, content=content)

    def upload_metadata(self, _: Sequence[Record]) -> None:
        metadata_path = self._run_dir / "files" / "swanlab-metadata.json"
        with safe.block(message=f"Failed to upload metadata, skipping; file kept at {metadata_path}"):
            with open(metadata_path, "r", encoding="utf-8") as f:
                content = json.load(f)
            if isinstance(content, dict):
                upload_metadata(self._username, self._project, self._cuid, content=content)

    def upload_requirements(self, _: Sequence[Record]) -> None:
        requirements_path = self._run_dir / "files" / "requirements.txt"
        with safe.block(message=f"Failed to upload requirements, skipping; file kept at {requirements_path}"):
            with open(requirements_path, "r", encoding="utf-8") as f:
                content = f.read()
            if len(content) > 0:
                upload_requirements(self._username, self._project, self._cuid, content=content)

    def upload_conda(self, _: Sequence[Record]) -> None:
        conda_path = self._run_dir / "files" / "conda.yaml"
        with safe.block(message=f"Failed to upload conda, skipping; file kept at {conda_path}"):
            with open(conda_path, "r", encoding="utf-8") as f:
                content = f.read()
            if len(content) > 0:
                upload_conda(self._username, self._project, self._cuid, content=content)

    def close(self) -> None:
        """关闭 sender。"""
        return


__all__ = [
    "HttpRecordSender",
]
