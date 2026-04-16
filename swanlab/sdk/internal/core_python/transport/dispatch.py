"""
@author: caddiesnew
@file: dispatch.py
@time: 2026/4/16
@description: Record 按类型聚合分发 + 失败回滚
"""

import threading
from typing import Callable, List, Optional

from swanlab.proto.swanlab.record.v1.record_pb2 import Record
from swanlab.sdk.internal.core_python.transport.helper import generate_chunks, group_records_by_type
from swanlab.sdk.internal.core_python.transport.sender import create_record_sender
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

    def __init__(
        self,
        cond: threading.Condition,
        buffer: List[Record],
        upload_callback: Optional[Callable[[int], None]] = None,
    ):
        self._cond = cond
        self._buffer = buffer
        self._upload_callback = upload_callback

    def __call__(self, records: List[Record]) -> None:
        """聚合 + 分发入口。"""
        grouped = group_records_by_type(records)
        for kind, typed_records in grouped.items():
            handler = getattr(self, f"_handle_{kind}", None)
            if handler:
                handler(typed_records)
            else:
                console.warning(f"No handler for record kind={kind!r}, skipping {len(typed_records)} records.")

    # ── 通用 record 上传函数 + 失败回滚 ──
    def _upload_typed(self, record_type: str, records: List[Record]) -> None:
        """
        将 records 按类型处理：分片 → sender 上传 → 回调。
        上传失败时回滚到 buffer 头部。
        sender 由模块级单例管理，无需在此关闭。
        """
        sender = create_record_sender()
        uploaded_count = 0
        try:
            for chunk, chunk_len in generate_chunks(records, _PER_REQUEST_LEN):
                sender.upload(record_type, chunk)
                uploaded_count += chunk_len
                if self._upload_callback:
                    self._upload_callback(chunk_len)
        except Exception:
            console.trace(f"upload error for record_type={record_type!r}")
            with self._cond:
                # 只回滚尚未成功上传的记录，避免重复上传
                self._buffer[:0] = records[uploaded_count:]

    # ── 类型 handler（按 _handle_{WhichOneof 返回值} 约定） ──

    def _handle_run(self, records: List[Record]) -> None:
        self._upload_typed("run", records)

    def _handle_finish(self, records: List[Record]) -> None:
        self._upload_typed("finish", records)

    def _handle_column(self, records: List[Record]) -> None:
        self._upload_typed("column", records)

    def _handle_metric(self, records: List[Record]) -> None:
        self._upload_typed("metric", records)

    def _handle_config(self, records: List[Record]) -> None:
        self._upload_typed("config", records)

    def _handle_console(self, records: List[Record]) -> None:
        self._upload_typed("console", records)

    def _handle_metadata(self, records: List[Record]) -> None:
        self._upload_typed("metadata", records)

    def _handle_requirements(self, records: List[Record]) -> None:
        self._upload_typed("requirements", records)

    def _handle_conda(self, records: List[Record]) -> None:
        self._upload_typed("conda", records)


__all__ = [
    "Dispatch",
]
