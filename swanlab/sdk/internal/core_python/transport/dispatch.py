"""
@author: caddiesnew
@file: dispatch.py
@time: 2026/4/16
@description: Record 按类型聚合分发 + 失败回滚
"""

import threading
import time
from typing import Callable, List, Optional, Sequence

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.helper import generate_chunks, group_records_by_type
from swanlab.sdk.internal.core_python.transport.sender import HttpRecordSender
from swanlab.sdk.internal.pkg import console, safe

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
    _RECORD_TYPES = frozenset(
        ("run", "finish", "column", "metric", "config", "console", "metadata", "requirements", "conda")
    )

    def __init__(
        self,
        cond: threading.Condition,
        buffer: List[Record],
        upload_callback: Optional[Callable[[int], None]] = None,
        max_retries: int = 3,
        initial_backoff: float = 0.5,
        sender: Optional[HttpRecordSender] = None,
    ):
        self._cond = cond
        self._buffer = buffer
        self._upload_callback = upload_callback
        self._max_retries = max_retries
        self._initial_backoff = initial_backoff
        self._sender: Optional[HttpRecordSender] = sender

    def __call__(self, records: List[Record]) -> None:
        """聚合 + 分发入口。"""
        grouped = group_records_by_type(records)
        for kind, typed_records in grouped.items():
            self._handle_record_by_type(kind, typed_records)

    # ── chunk 上传 ──
    @safe.decorator(level="error", message="record chunk upload failed")
    def _upload_chunk(self, record_type: str, chunk: Sequence[Record]) -> Optional[bool]:
        """上传单个 chunk。成功正常返回，失败时 safe.decorator 捕获异常并 console.trace。"""
        if self._sender is None:
            raise RuntimeError("sender not set")
        self._sender.upload(record_type, chunk)
        return True

    # ── 通用 record 上传函数 + 指数退避重试 + 失败回滚 ──
    def _upload_typed(self, record_type: str, records: List[Record]) -> None:
        """
        将 records 按类型处理：分片 → sender 上传 → 回调。
        上传失败时按指数退避重试，重试耗尽后回滚到 buffer 头部。
        """
        uploaded_count = 0
        for chunk, chunk_len in generate_chunks(records, _PER_REQUEST_LEN):
            success = False
            for attempt in range(self._max_retries):
                result = self._upload_chunk(record_type, chunk)
                if result is True:
                    success = True
                    break
                # 指数退避：最后一次失败不等待
                if attempt < self._max_retries - 1:
                    delay = self._initial_backoff * (2**attempt)
                    time.sleep(delay)

            if success:
                uploaded_count += chunk_len
                if self._upload_callback:
                    self._upload_callback(chunk_len)
            else:
                # 重试耗尽 → 回滚未上传记录到 buffer，然后 return
                with self._cond:
                    self._buffer[:0] = records[uploaded_count:]
                return

    def _handle_record_by_type(self, kind: str, records: List[Record]) -> None:
        """按 kind 路由到 _upload_typed，与 HttpRecordSender.upload_{kind} 对齐。"""
        if kind in self._RECORD_TYPES:
            self._upload_typed(kind, records)
        else:
            console.warning(f"No handler for record kind={kind!r}, skipping {len(records)} records.")


__all__ = [
    "Dispatch",
]
