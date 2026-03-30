"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传入口与 Core sidecar transport
"""

import os
import time
from dataclasses import dataclass
from typing import Callable, Optional, Protocol, Sequence

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import console

from .helper import generate_chunks, group_records_by_type

_TRUE_VALUES = {"true", "1", "yes", "on"}
_FALSE_VALUES = {"false", "0", "no", "off"}
DEFAULT_CORE_TIMEOUT = 10.0


def _parse_bool(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    lowered = raw.strip().lower()
    if lowered in _TRUE_VALUES:
        return True
    if lowered in _FALSE_VALUES:
        return False
    console.warning(f"Invalid {name} value: {raw!r}. Using default {default}.")
    return default


def _parse_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    try:
        value = int(raw)
    except ValueError:
        console.warning(f"Invalid {name} value: {raw!r}. Using default {default}.")
        return default
    if value < 0:
        console.warning(f"{name} must be >= 0, got {value}. Using default {default}.")
        return default
    return value


def _parse_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    try:
        value = float(raw)
    except ValueError:
        console.warning(f"Invalid {name} value: {raw!r}. Using default {default}.")
        return default
    if value <= 0:
        console.warning(f"{name} must be > 0, got {value}. Using default {default}.")
        return default
    return value


@dataclass(frozen=True)
class CoreTransportConfig:
    """
    Python -> Core sidecar 传输配置。

    这里只负责环境变量解析；真实的 Core RPC 实现后续再接入。
    """

    enabled: bool
    address: Optional[str]
    timeout: float

    @classmethod
    def from_env(cls) -> "CoreTransportConfig":
        address = os.getenv("SWANLAB_CORE_ADDRESS") or None
        host = os.getenv("SWANLAB_CORE_HOST") or None
        port = _parse_int("SWANLAB_CORE_PORT", 0)

        if address is None and host and port > 0:
            address = f"{host}:{port}"

        enabled = _parse_bool("SWANLAB_CORE_ENABLED", default=True)
        timeout = _parse_float("SWANLAB_CORE_TIMEOUT", DEFAULT_CORE_TIMEOUT)
        return cls(enabled=enabled, address=address, timeout=timeout)


class RecordTransport(Protocol):
    """Record 上传 transport 协议，屏蔽具体 HTTP/REST 实现。"""

    def upload_record_group(self, record_type: str, records: Sequence[Record]) -> None:
        """按 record_type 批量上传一组 protobuf Record。"""
        ...

    def close(self) -> None:
        """关闭 transport，释放底层资源。"""
        ...


class HttpRecordTransport:
    """
    基于现有 HTTP client 的上传 transport 骨架。

    当前阶段只负责按 record_type 分发，为后续 REST body 构造预留扩展点。
    """

    _announced = False

    def __init__(self, config: CoreTransportConfig):
        self._config = config
        self._announce_once()

    def _announce_once(self) -> None:
        if HttpRecordTransport._announced:
            return
        if not self._config.enabled:
            console.debug("Core record transport is disabled; uploader transport is running in no-op mode.")
        else:
            console.debug("Core record transport is enabled and using the HTTP uploader skeleton.")
        HttpRecordTransport._announced = True

    def upload_record_group(self, record_type: str, records: Sequence[Record]) -> None:
        if len(records) == 0:
            return

        handler = getattr(self, f"_upload_{record_type}_group", None)
        if handler is None:
            self._upload_unknown_group(record_type, records)
            return
        handler(records)

    def _upload_metric_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("metric", records)

    def _upload_config_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("config", records)

    def _upload_run_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("run", records)

    def _upload_finish_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("finish", records)

    def _upload_column_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("column", records)

    def _upload_console_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("console", records)

    def _upload_metadata_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("metadata", records)

    def _upload_requirements_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("requirements", records)

    def _upload_conda_group(self, records: Sequence[Record]) -> None:
        self._upload_group_skeleton("conda", records)

    def _upload_unknown_group(self, record_type: str, records: Sequence[Record]) -> None:
        self._upload_group_skeleton(record_type, records)

    def _upload_group_skeleton(self, record_type: str, records: Sequence[Record]) -> None:
        del records
        if not self._config.enabled:
            return
        console.debug(f"HTTP upload skeleton is ready for record_type={record_type!r}, but request mapping is pending.")

    def close(self) -> None:
        pass


class NoopRecordTransport:
    """禁用上传时使用的空实现 transport。"""

    def __init__(self, config: CoreTransportConfig):
        self._config = config

    def upload_record_group(self, record_type: str, records: Sequence[Record]) -> None:
        del record_type, records
        if self._config.enabled:
            console.debug("NoopRecordTransport received upload request while enabled; request was ignored.")

    def close(self) -> None:
        pass


def create_record_transport(config: Optional[CoreTransportConfig] = None) -> RecordTransport:
    """
    创建 Record transport。

    当前阶段在启用时返回 HTTP transport 骨架，禁用时返回 no-op transport。
    """

    resolved = config or CoreTransportConfig.from_env()
    if not resolved.enabled:
        return NoopRecordTransport(resolved)
    return HttpRecordTransport(resolved)


def trace_records(
    records: Optional[Sequence[Record]] = None,
    per_request_len: int = 10_000,
    upload_callback: Optional[Callable[[int], None]] = None,
    transport: Optional[RecordTransport] = None,
) -> None:
    """
    分片上传 protobuf Record。

    线程层负责缓冲、聚合与重试；此层只做 Record 分片、按 record_type 分组和 transport 调用。
    """
    if records is None or len(records) == 0:
        return

    is_split_mode = per_request_len != -1 and len(records) > per_request_len
    owns_transport = transport is None
    active_transport = transport or create_record_transport()

    try:
        for chunk, chunk_len in generate_chunks(records, per_request_len):
            for record_type, grouped_records in group_records_by_type(chunk).items():
                active_transport.upload_record_group(record_type, grouped_records)
            if upload_callback:
                upload_callback(chunk_len)
            if is_split_mode:
                time.sleep(1)
    finally:
        if owns_transport:
            active_transport.close()


def upload_records(
    records: Sequence[Record],
    upload_callback: Optional[Callable[[int], None]] = None,
    per_request_len: int = 1000,
) -> None:
    """上传一组 protobuf Record，不在此层做任何 JSON 化。"""
    if len(records) == 0:
        return console.debug("No records to upload.")
    trace_records(records, per_request_len=per_request_len, upload_callback=upload_callback)


__all__ = [
    "CoreTransportConfig",
    "HttpRecordTransport",
    "NoopRecordTransport",
    "RecordTransport",
    "create_record_transport",
    "trace_records",
    "upload_records",
]
