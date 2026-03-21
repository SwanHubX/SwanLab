"""
@author: caddiesnew
@file: upload.py
@time: 2026/3/19
@description: Record 上传入口
"""

from typing import Callable, Optional, Sequence

from swanlab.sdk.internal.pkg import console

from .batch import RecordLike, load_record, trace_records


def upload_records(
    records: Sequence[RecordLike],
    upload_callback: Optional[Callable[[int], None]] = None,
    per_request_len: int = 1000,
) -> None:
    """上传一组 protobuf Record，不在此层做任何 JSON 化。"""
    if len(records) == 0:
        return console.debug("No records to upload.")
    trace_records(records, per_request_len=per_request_len, upload_callback=upload_callback)


__all__ = ["RecordLike", "load_record", "trace_records", "upload_records"]
