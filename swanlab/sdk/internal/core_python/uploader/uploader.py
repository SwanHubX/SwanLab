"""
@author: caddiesnew
@file: uploader.py
@time: 2026/3/29
@description: Timer-driven protobuf batch uploader
"""

import threading
from queue import Empty

from swanlab.sdk.internal.bus.events import FlushPayload, MetricsUploadEvent
from swanlab.sdk.internal.pkg.timer import Timer

from swanlab.sdk.internal.core_python.uploader.helper import RecordQueue
from swanlab.sdk.internal.core_python.uploader.sender import HttpRecordSender, group_records_by_type


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
        self._queue = RecordQueue()
        self._pending: FlushPayload = []
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

    def enqueue(self, records: FlushPayload) -> None:
        """Append a protobuf record batch into the uploader queue."""
        if self._closed:
            raise RuntimeError("HttpBatchUploader is closed.")
        if len(records) == 0:
            return
        self._queue.put_all(records)

    def _drain_queue(self) -> None:
        while True:
            try:
                self._pending.append(self._queue.get_nowait())
            except Empty:
                return

    def flush(self) -> None:
        """Send the current queued records grouped by record_type."""
        with self._flush_lock:
            self._drain_queue()
            if len(self._pending) == 0:
                return

            buckets = group_records_by_type(self._pending)
            self._sender.send(MetricsUploadEvent(buckets=buckets))
            self._pending = []

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
