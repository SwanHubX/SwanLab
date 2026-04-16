"""
@author: caddiesnew
@file: sender.py
@time: 2026/3/19
@description: Record 上传 transport 抽象（HTTP 传输层）
"""

from typing import Sequence

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.pkg import console


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


__all__ = [
    "HttpRecordTransport",
    "create_record_transport",
]
