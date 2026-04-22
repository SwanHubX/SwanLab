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
from swanlab.sdk.internal.pkg import console


class HttpRecordSender:
    """HTTP 传输层上传 sender，请求映射待实现。"""

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

    def upload_config(self, records: Sequence[Record]) -> None:
        with open(self._run_dir / "files" / "config.yaml", "r", encoding="utf-8") as f:
            content = yaml.safe_load(f)
        upload_config(self._username, self._project, self._cuid, content=content)

    def upload_console(self, records: Sequence[Record]) -> None:
        console.debug("HTTP upload skeleton: upload_console (request mapping pending).")

    def upload_metadata(self, records: Sequence[Record]) -> None:
        with open(self._run_dir / "files" / "swanlab-metadata.json", "r", encoding="utf-8") as f:
            content = json.load(f)
        upload_metadata(self._username, self._project, self._cuid, content=content)

    def upload_requirements(self, records: Sequence[Record]) -> None:
        with open(self._run_dir / "files" / "requirements.txt", "r", encoding="utf-8") as f:
            content = f.read()
        upload_requirements(self._username, self._project, self._cuid, content=content)

    def upload_conda(self, records: Sequence[Record]) -> None:
        with open(self._run_dir / "files" / "conda.yaml", "r", encoding="utf-8") as f:
            content = f.read()
        upload_conda(self._username, self._project, self._cuid, content=content)

    def close(self) -> None:
        """关闭 sender。"""
        return


__all__ = [
    "HttpRecordSender",
]
