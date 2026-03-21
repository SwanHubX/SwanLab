"""
@author: caddiesnew
@file: upload.py
@time: 2026/3/19
@description: Record 上传入口与 Core sidecar transport
"""

import os
import time
from dataclasses import dataclass
from typing import Callable, Optional, Protocol, Sequence

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import console

from .batch import RecordLike, _generate_chunks, load_record

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

        enabled = _parse_bool("SWANLAB_CORE_ENABLED", default=bool(address))
        timeout = _parse_float("SWANLAB_CORE_TIMEOUT", DEFAULT_CORE_TIMEOUT)
        return cls(enabled=enabled, address=address, timeout=timeout)


class RecordTransport(Protocol):
    """Record 上传 transport 协议，屏蔽具体 RPC 实现。"""

    def upsert_record(self, record: Record) -> None:
        """上传单条 protobuf Record。"""
        ...

    def close(self) -> None:
        """关闭 transport，释放底层资源。"""
        ...


class NoopRecordTransport:
    """
    Core sidecar 尚未落地前的占位 transport。

    线程层仍可完成缓冲、聚合和重试，但此处不做真实网络发送。
    """

    _announced = False

    def __init__(self, config: CoreTransportConfig):
        self._config = config
        self._announce_once()

    def _announce_once(self) -> None:
        if NoopRecordTransport._announced:
            return
        if not self._config.enabled:
            console.debug("Core record transport is disabled; uploader transport is running in no-op mode.")
        elif self._config.address:
            console.debug(
                f"Core record transport is configured for {self._config.address}, "
                "but the sidecar transport is not implemented yet."
            )
        else:
            console.debug("Core record transport is enabled, but no Core address is configured yet.")
        NoopRecordTransport._announced = True

    def upsert_record(self, record: Record) -> None:
        del record
        pass

    def close(self) -> None:
        pass


def create_record_transport(config: Optional[CoreTransportConfig] = None) -> RecordTransport:
    """
    创建 Record transport。

    当前阶段固定返回 no-op placeholder，待 Core sidecar 落地后在此替换成真实实现。
    """

    resolved = config or CoreTransportConfig.from_env()
    return NoopRecordTransport(resolved)


def trace_records(
    records: Optional[Sequence[RecordLike]] = None,
    per_request_len: int = 1000,
    upload_callback: Optional[Callable[[int], None]] = None,
    transport: Optional[RecordTransport] = None,
) -> None:
    """
    分片上传 protobuf Record。

    线程层负责缓冲、聚合与重试；此层只做 Record 反序列化、分片和 transport 调用。
    """
    if records is None or len(records) == 0:
        return

    is_split_mode = per_request_len != -1 and len(records) > per_request_len
    owns_transport = transport is None
    active_transport = transport or create_record_transport()

    try:
        for chunk, chunk_len in _generate_chunks(records, per_request_len):
            for record_like in chunk:
                active_transport.upsert_record(load_record(record_like))
            if upload_callback:
                upload_callback(chunk_len)
            if is_split_mode:
                time.sleep(1)
    finally:
        if owns_transport:
            active_transport.close()


def upload_records(
    records: Sequence[RecordLike],
    upload_callback: Optional[Callable[[int], None]] = None,
    per_request_len: int = 1000,
) -> None:
    """上传一组 protobuf Record，不在此层做任何 JSON 化。"""
    if len(records) == 0:
        return console.debug("No records to upload.")
    trace_records(records, per_request_len=per_request_len, upload_callback=upload_callback)


__all__ = [
    "CoreTransportConfig",
    "NoopRecordTransport",
    "RecordLike",
    "RecordTransport",
    "create_record_transport",
    "load_record",
    "trace_records",
    "upload_records",
]
