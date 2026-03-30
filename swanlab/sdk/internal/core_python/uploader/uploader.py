"""
@author: caddiesnew
@file: uploader.py
@time: 2026/3/29
@description: Timer-driven protobuf batch uploader
"""

import threading
from queue import Empty

from swanlab.sdk.internal.bus.events import FlushPayload, MetricsUploadEvent
from swanlab.sdk.internal.core_python.uploader.helper import RecordQueue
from swanlab.sdk.internal.core_python.uploader.sender import HttpRecordSender, group_records_by_type
from swanlab.sdk.internal.pkg import console
from swanlab.sdk.internal.pkg.timer import Timer


class HttpBatchUploader:
    """Buffer protobuf records in memory and flush them to an HTTP sender periodically."""

    TIMER_NAME = "SwanLab·Uploader"
    DEFAULT_MAX_PENDING = 100_000

    def __init__(
        self,
        sender: HttpRecordSender,
        upload_interval: float = 1.0,
        max_pending: int = DEFAULT_MAX_PENDING,
        auto_start: bool = True,
    ) -> None:
        if max_pending <= 0:
            raise ValueError(f"max_pending must be greater than 0, got {max_pending}.")

        self._sender = sender
        self._queue = RecordQueue()
        self._pending: FlushPayload = []
        self._flush_lock = threading.Lock()
        self._state_lock = threading.Lock()
        self._closed = False
        self._max_pending = max_pending
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
        with self._state_lock:
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

    def _trim_pending(self) -> None:
        overflow = len(self._pending) - self._max_pending
        if overflow <= 0:
            return
        console.warning(
            f"Uploader pending buffer exceeded limit ({self._max_pending}), dropping {overflow} oldest records."
        )
        self._pending = self._pending[overflow:]

    def flush(self) -> None:
        """Send the current queued records grouped by record_type."""
        with self._flush_lock:
            self._drain_queue()
            if len(self._pending) == 0:
                return
            self._trim_pending()

            buckets = group_records_by_type(self._pending)
            self._sender.send(MetricsUploadEvent(buckets=buckets))
            self._pending = []

    def close(self) -> None:
        """Stop periodic flush and perform a final synchronous flush."""
        with self._state_lock:
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
