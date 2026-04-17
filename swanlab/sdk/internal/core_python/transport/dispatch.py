"""
@author: caddiesnew
@file: dispatch.py
@time: 2026/4/16
@description: Record 按类型聚合分发 + 失败回滚
"""

import time
from typing import Callable, List, Optional, Tuple

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.helper import generate_chunks, group_records_by_type
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.pkg import console

# 单次请求最大 record 数
_MAX_RECORDS_PER_REQUEST = 10_000


class Dispatch:
    """
    按 record_type 聚合 records 并分发给对应 handler。

    遍历 pending records，用 WhichOneof("record_type") 取 record_type 并按类型分组，
    再逐个交给 _upload_record_type() 上传。
    上传失败时将当前类型剩余 records 与后续未处理类型一起返回给上层保留。
    """

    # ── record 类型分发 ──
    _RECORD_TYPES = frozenset(f.name for f in Record.DESCRIPTOR.oneofs_by_name["record_type"].fields)

    def __init__(
        self,
        upload_callback: Optional[Callable[[int], None]] = None,
        max_retries: int = 3,
        initial_backoff: float = 0.5,
        sender: Optional[HttpRecordSender] = None,
    ):
        self._upload_callback = upload_callback
        self._max_retries = max_retries
        self._initial_backoff = initial_backoff
        self._sender: Optional[HttpRecordSender] = sender

    def __call__(self, records: List[Record]) -> Tuple[bool, List[Record]]:
        """按 record_type 聚合并上传，返回 (是否成功, 失败及未处理的 records)。"""
        records_by_type = group_records_by_type(records)
        record_types = list(records_by_type.keys())

        for index, record_type in enumerate(record_types):
            success, remaining_records = self._upload_record_type(record_type, records_by_type[record_type])
            if not success:
                unprocessed_records = [
                    record
                    for next_record_type in record_types[index + 1 :]
                    for record in records_by_type[next_record_type]
                ]
                return False, remaining_records + unprocessed_records
        return True, []

    def _upload_record_type(self, record_type: str, records: List[Record]) -> Tuple[bool, List[Record]]:
        """上传单个 record_type 分组，返回 (是否成功, 当前分组剩余的 records)。"""
        if record_type not in self._RECORD_TYPES:
            console.warning(f"No handler for record kind={record_type!r}, skipping {len(records)} records.")
            return True, []

        uploaded_record_count = 0
        for chunk, chunk_size in generate_chunks(records, _MAX_RECORDS_PER_REQUEST):
            success = False
            for attempt in range(self._max_retries):
                try:
                    if self._sender is None:
                        raise RuntimeError("sender not set")
                    self._sender.upload(record_type, chunk)
                    success = True
                    break
                except Exception:
                    console.trace(f"record chunk upload failed, retrying ({attempt + 1}) times..", level_name="error")
                    if attempt < self._max_retries - 1:
                        time.sleep(self._initial_backoff * (2**attempt))

            if success:
                uploaded_record_count += chunk_size
                if self._upload_callback:
                    self._upload_callback(chunk_size)
            else:
                return False, records[uploaded_record_count:]
        return True, []


__all__ = [
    "Dispatch",
]
