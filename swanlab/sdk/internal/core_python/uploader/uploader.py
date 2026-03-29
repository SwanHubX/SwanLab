"""
@author: caddiesnew
@file: uploader.py
@time: 2026/3/29
@description: Timer-driven protobuf batch uploader
"""

import threading

from swanlab.sdk.internal.bus.events import MetricsUploadEvent
from swanlab.sdk.internal.pkg.timer import Timer
from swanlab.proto.swanlab.record.v1.record_pb2 import Record

from typing import List

from .http import HttpRecordSender, group_records_by_type


class HttpBatchUploader:
    """Buffer protobuf records in memory and flush them to an HTTP sender periodically."""

    TIMER_NAME = "SwanLab·Uploader"

    def __init__(
        self,
        sender: HttpRecordSender,
        upload_interval: float = 1.0,
        auto_start: bool = True,
    ) -> None:
        self._sender = sender
        self._buffer: List[Record] = []
        self._buffer_lock = threading.Lock()
        self._flush_lock = threading.Lock()
        self._closed = False
        self._timer = Timer(
            self.flush,
            interval=upload_interval,
            immediate=False,
            name=self.TIMER_NAME,
        )

        if auto_start:
            self._timer.start()

    def enqueue(self, records: List[Record]) -> None:
        """Append a protobuf record batch into the in-memory upload buffer."""
        if self._closed:
            raise RuntimeError("HttpBatchUploader is closed.")
        if len(records) == 0:
            return
        with self._buffer_lock:
            self._buffer.extend(records)

    def flush(self) -> None:
        """Send the current buffered records grouped by record_type."""
        with self._flush_lock:
            with self._buffer_lock:
                if len(self._buffer) == 0:
                    return
                pending = list(self._buffer)

            buckets = group_records_by_type(pending)
            self._sender.send(MetricsUploadEvent(buckets=buckets))

            with self._buffer_lock:
                del self._buffer[: len(pending)]

    def close(self) -> None:
        """Stop periodic flush and perform a final synchronous flush."""
        if self._closed:
            return

        self._closed = True
        self._timer.cancel()
        self._timer.join(timeout=10)
        try:
            self.flush()
        finally:
            self._sender.close()


__all__ = [
    "HttpBatchUploader",
]
