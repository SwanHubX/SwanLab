"""
@author: caddiesnew
@file: batch.py
@time: 2026/3/19
@description: 基于 gRPC 的 Record 上传传输层
"""

import time
from typing import Callable, Iterator, List, Optional, Sequence, Tuple, Union
from urllib.parse import urlparse
from weakref import WeakKeyDictionary

import grpc

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.proto.swanlab.record.v1.record_pb2_grpc import RecordServiceStub
from swanlab.sdk.internal.core_python.client import Client
from swanlab.sdk.internal.core_python.client import _get_client as get_client

RecordLike = Union[Record, bytes]

_CHANNELS: "WeakKeyDictionary[Client, grpc.Channel]" = WeakKeyDictionary()
_STUBS: "WeakKeyDictionary[Client, RecordServiceStub]" = WeakKeyDictionary()
DEFAULT_RPC_TIMEOUT = 10.0


def load_record(record_like: RecordLike) -> Record:
    """统一接收 protobuf Record 或其序列化 bytes。"""
    if isinstance(record_like, Record):
        return record_like
    if isinstance(record_like, bytes):
        record = Record()
        record.ParseFromString(record_like)
        return record
    raise TypeError(f"Unsupported record type: {type(record_like).__name__}")


def _generate_chunks(records: Sequence[RecordLike], per_request_len: int) -> Iterator[Tuple[Sequence[RecordLike], int]]:
    """
    按固定大小对 Record 序列做切片。
    yield: (分片数据, 分片长度)
    """
    if per_request_len == -1:
        yield records, len(records)
        return

    total = len(records)
    if total <= per_request_len:
        yield records, total
        return

    for index in range(0, total, per_request_len):
        chunk = records[index : index + per_request_len]
        yield chunk, len(chunk)


def _build_rpc_target(client: Client) -> Tuple[str, bool]:
    parsed = urlparse(client._base_url)
    hostname = parsed.hostname
    if hostname is None:
        raise RuntimeError(f"Invalid API host for gRPC transport: {client._base_url}")

    secure = parsed.scheme != "http"
    port = parsed.port or (443 if secure else 80)
    return f"{hostname}:{port}", secure


def _get_channel(client: Client) -> grpc.Channel:
    channel = _CHANNELS.get(client)
    if channel is not None:
        return channel

    target, secure = _build_rpc_target(client)
    options = [
        ("grpc.keepalive_time_ms", 30_000),
        ("grpc.keepalive_timeout_ms", 10_000),
    ]
    if secure:
        channel = grpc.secure_channel(target, grpc.ssl_channel_credentials(), options=options)
    else:
        channel = grpc.insecure_channel(target, options=options)
    _CHANNELS[client] = channel
    return channel


def _get_record_service(client: Client) -> RecordServiceStub:
    stub = _STUBS.get(client)
    if stub is not None:
        return stub

    stub = RecordServiceStub(_get_channel(client))
    _STUBS[client] = stub
    return stub


def _build_rpc_metadata(client: Client) -> List[Tuple[str, str]]:
    # 复用现有 Client 的鉴权刷新逻辑，确保 sid 在 gRPC 调用前仍然有效。
    client._before_request()
    sid = client._session.cookies.get("sid")
    if sid is None:
        raise RuntimeError("SwanLab client is not authenticated.")
    return [
        ("cookie", f"sid={sid}"),
        ("x-swanlab-sdk-version", client._version),
    ]


def trace_records(
    records: Optional[Sequence[RecordLike]] = None,
    per_request_len: int = 1000,
    upload_callback: Optional[Callable] = None,
    timeout: float = DEFAULT_RPC_TIMEOUT,
):
    """
    分片调用 RecordService.UpsertRecord 上传 protobuf Record。
    :param records: 要上传的 Record 集合
    :param per_request_len: 每批上传的数量
    :param upload_callback: 上传进度回调函数
    :param timeout: 单次 RPC 超时时间（秒）
    """
    if records is None or len(records) == 0:
        return

    is_split_mode = per_request_len != -1 and len(records) > per_request_len
    client = get_client()
    stub = _get_record_service(client)

    for chunk, chunk_len in _generate_chunks(records, per_request_len):
        metadata = _build_rpc_metadata(client)
        for record_like in chunk:
            stub.UpsertRecord(load_record(record_like), timeout=timeout, metadata=metadata, wait_for_ready=True)
        if upload_callback:
            upload_callback(chunk_len)
        if is_split_mode:
            time.sleep(1)


__all__ = [
    "RecordLike",
    "load_record",
    "trace_records",
]
