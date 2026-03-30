"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传入口与 Core sidecar transport
"""

import time
from typing import Callable, Optional, Sequence

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import console

from .helper import generate_chunks, group_records_by_type


class HttpRecordTransport:
    """PLACEHOLDER: 基于现有 HTTP client 的上传 transport 。"""

    def upload_record_group(self, record_type: str, records: Sequence[Record]) -> None:
        if len(records) == 0:
            return
        console.debug(f"HTTP upload skeleton is ready for record_type={record_type!r}, but request mapping is pending.")

    def close(self) -> None:
        pass


def create_record_transport() -> HttpRecordTransport:
    """创建 Record transport。"""
    return HttpRecordTransport()


def trace_records(
    records: Optional[Sequence[Record]] = None,
    per_request_len: int = 10_000,
    upload_callback: Optional[Callable[[int], None]] = None,
    transport: Optional[HttpRecordTransport] = None,
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
    "HttpRecordTransport",
    "create_record_transport",
    "trace_records",
    "upload_records",
]
