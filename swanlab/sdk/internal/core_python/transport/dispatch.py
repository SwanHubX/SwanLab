"""
@author: caddiesnew
@file: dispatch.py
@time: 2026/4/16
@description: Record 按类型聚合分发 + 失败回滚
"""

import threading
import time
from typing import Callable, List, Optional, Tuple

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.buffer import RecordBuffer
from swanlab.sdk.internal.core_python.transport.helper import generate_chunks, group_records_by_type
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.pkg import console

# 单次请求最大 record 数
_PER_REQUEST_LEN = 10_000


class Dispatch:
    """
    按 record_type 聚合 records 并分发给对应 handler。

    遍历 pending records，用 WhichOneof("record_type") 取 kind 并按 kind 分组，
    再通过 _handle_{kind} 约定路由到对应 handler。
    上传失败时通过 Condition 锁回滚到 buffer 头部。
    """

    # ── record 类型分发 ──
    _RECORD_TYPES = frozenset(f.name for f in Record.DESCRIPTOR.oneofs_by_name["record_type"].fields)

    def __init__(
        self,
        cond: threading.Condition,
        buf: RecordBuffer,
        upload_callback: Optional[Callable[[int], None]] = None,
        max_retries: int = 3,
        initial_backoff: float = 0.5,
        sender: Optional[HttpRecordSender] = None,
    ):
        self._cond = cond
        self._buf = buf
        self._upload_callback = upload_callback
        self._max_retries = max_retries
        self._initial_backoff = initial_backoff
        self._sender: Optional[HttpRecordSender] = sender

    def __call__(self, records: List[Record]) -> bool:
        """聚合 + 分发入口。返回 True 表示全部成功，False 表示有回滚。"""
        grouped = group_records_by_type(records)
        keys = list(grouped.keys())
        for i, kind in enumerate(keys):
            success, failed = self._handle_record_by_type(kind, grouped[kind])
            if not success:
                # 当前组失败部分 + 后续所有未处理组统一回滚，保持顺序
                to_rollback = failed + [record for key in keys[i + 1 :] for record in grouped[key]]
                with self._cond:
                    self._buf.prepend(to_rollback)
                return False
        return True

    def _handle_record_by_type(self, kind: str, records: List[Record]) -> Tuple[bool, List[Record]]:
        """按 kind 上传 record，指数退避重试。返回 (是否成功, 失败的 records)。"""
        if kind not in self._RECORD_TYPES:
            console.warning(f"No handler for record kind={kind!r}, skipping {len(records)} records.")
            return True, []

        uploaded_count = 0
        for chunk, chunk_len in generate_chunks(records, _PER_REQUEST_LEN):
            success = False
            for attempt in range(self._max_retries):
                try:
                    if self._sender is None:
                        raise RuntimeError("sender not set")
                    self._sender.upload(kind, chunk)
                    success = True
                    break
                except Exception:
                    console.trace("record chunk upload failed", level_name="error")
                    if attempt < self._max_retries - 1:
                        time.sleep(self._initial_backoff * (2**attempt))

            if success:
                uploaded_count += chunk_len
                if self._upload_callback:
                    self._upload_callback(chunk_len)
            else:
                return False, records[uploaded_count:]
        return True, []


__all__ = [
    "Dispatch",
]
